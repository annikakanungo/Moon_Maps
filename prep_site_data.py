#!/usr/bin/env python3
"""
Build terrain-aware Sun/Earth visibility data for lunar south pole sites.

For each site this writes data/<id>.json containing
  - horizon: terrain horizon elevation (deg) for every azimuth 0..359 (from LRO LOLA DEM)
  - sun / earth: hourly azimuth and elevation (deg) from NASA JPL Horizons
and prints % sunlit, % Earth in view and % both, using the real horizon.

Setup (once):
  pip install numpy rasterio pyproj
  Download the 80 m south polar LOLA DEM from the PDS Geosciences Node (LOLA gridded
  data records, file LDEM_80S_80MPP_ADJ.TIF) and put it at data/LDEM_80S_80MPP_ADJ.TIF
  The DEM must be in metres, polar stereographic, covering 80S to 90S.

Run:
  python prep_site_data.py --start 2026-11-14 --days 30
"""
import argparse, json, math, os, re, urllib.parse, urllib.request
import numpy as np
import rasterio
from pyproj import CRS, Transformer
from rasterio.windows import Window

R = 1737400.0  # lunar reference radius in metres (matches the LOLA datum)
SUN_RADIUS_DEG = 0.267

SITES = {  # id: (name, latitude, east longitude)
    "shackleton": ("Shackleton Crater (south pole)", -89.9, 0.0),
    "malapert": ("Malapert Massif", -85.99, 0.12),
    "amundsen_rim": ("Amundsen Rim (site AR01)", -84.467, 66.97),
    "nobile_rim_2": ("Nobile Rim 2 (study site)", -84.2016, 60.6999),
}
DEFAULT_PROJ4 = "+proj=stere +lat_0=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=1737400 +b=1737400 +units=m +no_defs"


def destination(lat, lon, az_deg, dist_m):
    """Point reached from (lat, lon) travelling dist_m along azimuth az_deg on a sphere."""
    p, l, a = map(math.radians, (lat, lon, az_deg))
    d = dist_m / R
    p2 = np.arcsin(np.sin(p) * np.cos(d) + np.cos(p) * np.sin(d) * np.cos(a))
    l2 = l + np.arctan2(np.sin(a) * np.sin(d) * np.cos(p), np.cos(d) - np.sin(p) * np.sin(p2))
    return np.degrees(p2), (np.degrees(l2) + 180) % 360 - 180


def horizon_mask(dem_path, lat, lon, max_km=30.0, mast_m=0.0, proj4=DEFAULT_PROJ4):
    """Highest terrain elevation angle (deg) seen from the site, for each azimuth 0..359."""
    geo = CRS.from_proj4("+proj=longlat +a=1737400 +b=1737400 +no_defs")
    with rasterio.open(dem_path) as src:
        crs = src.crs or CRS.from_proj4(proj4)
        tr = Transformer.from_crs(geo, crs, always_xy=True)
        res = abs(src.transform.a)
        pad = int(max_km * 1000 / res) + 5
        x0, y0 = tr.transform(lon, lat)
        c0, r0 = ~src.transform * (x0, y0)
        win = Window(int(c0) - pad, int(r0) - pad, 2 * pad + 1, 2 * pad + 1)
        z = src.read(1, window=win, boundless=True, fill_value=np.nan).astype("float64")
        wtf = src.window_transform(win)
        nodata = src.nodata
    if nodata is not None:
        z[z == nodata] = np.nan
    inv = ~wtf

    def height(la, lo):
        x, y = tr.transform(lo, la)
        c, r = inv * (x, y)
        c, r = np.asarray(c), np.asarray(r)
        c0_, r0_ = np.floor(c).astype(int), np.floor(r).astype(int)
        fc, fr = c - c0_, r - r0_
        ok = (c0_ >= 0) & (r0_ >= 0) & (c0_ < z.shape[1] - 1) & (r0_ < z.shape[0] - 1)
        c0_, r0_ = np.clip(c0_, 0, z.shape[1] - 2), np.clip(r0_, 0, z.shape[0] - 2)
        h = (z[r0_, c0_] * (1 - fc) * (1 - fr) + z[r0_, c0_ + 1] * fc * (1 - fr)
             + z[r0_ + 1, c0_] * (1 - fc) * fr + z[r0_ + 1, c0_ + 1] * fc * fr)
        return np.where(ok, h, np.nan)

    h0 = float(height(np.array([lat]), np.array([lon]))[0]) + mast_m
    dist = np.arange(res, max_km * 1000, res)
    delta = dist / R
    out = []
    for az in range(360):
        la, lo = destination(lat, lon, az, dist)
        h = height(la, lo)
        # exact spherical elevation angle, includes lunar curvature
        el = np.degrees(np.arctan2((R + h) * np.cos(delta) - (R + h0), (R + h) * np.sin(delta)))
        out.append(float(np.nanmax(el)) if np.isfinite(el).any() else 0.0)
    return out


def horizons(body, lat, lon, start, stop, step="1h"):
    """Hourly apparent azimuth/elevation of a body (10 = Sun, 399 = Earth) from a lunar surface point."""
    q = {
        "format": "text", "COMMAND": f"'{body}'", "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'OBSERVER'", "CENTER": "'coord@301'", "COORD_TYPE": "'GEODETIC'",
        "SITE_COORD": f"'{lon},{lat},0'", "START_TIME": f"'{start}'", "STOP_TIME": f"'{stop}'",
        "STEP_SIZE": f"'{step}'", "QUANTITIES": "'4'", "ANG_FORMAT": "'DEG'", "CSV_FORMAT": "'YES'",
    }
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(q)
    return parse_horizons(urllib.request.urlopen(url, timeout=120).read().decode())


def parse_horizons(text):
    m = re.search(r"\$\$SOE(.*?)\$\$EOE", text, re.S)
    if not m:
        raise RuntimeError("Horizons returned no ephemeris:\n" + text[:800])
    az, el = [], []
    for line in m.group(1).strip().splitlines():
        nums = []
        for part in line.split(",")[1:]:
            try:
                nums.append(float(part))
            except ValueError:
                pass  # flag columns such as '*' or 'C'
        az.append(nums[-2]); el.append(nums[-1])
    return {"az": az, "el": el}


def summarise(rec):
    hz = rec["horizon"]
    n = len(rec["sun"]["el"])
    sun = [rec["sun"]["el"][i] > hz[int(round(rec["sun"]["az"][i])) % 360] - SUN_RADIUS_DEG for i in range(n)]
    ear = [rec["earth"]["el"][i] > hz[int(round(rec["earth"]["az"][i])) % 360] for i in range(n)]
    pct = lambda v: 100 * sum(v) / n
    return pct(sun), pct(ear), pct([a and b for a, b in zip(sun, ear)])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dem", default="data/LDEM_80S_80MPP_ADJ.TIF")
    ap.add_argument("--start", required=True, help="UTC date, e.g. 2026-11-14")
    ap.add_argument("--days", type=int, default=30)
    ap.add_argument("--max-km", type=float, default=30.0)
    ap.add_argument("--mast", type=float, default=2.0, help="antenna/array height in metres")
    a = ap.parse_args()
    os.makedirs("data", exist_ok=True)
    stop = (np.datetime64(a.start) + np.timedelta64(a.days, "D")).astype(str)
    for sid, (name, lat, lon) in SITES.items():
        print(f"{name}: terrain horizon ...", flush=True)
        rec = {"id": sid, "name": name, "lat": lat, "lon": lon, "start": a.start, "days": a.days,
               "horizon": horizon_mask(a.dem, lat, lon, a.max_km, a.mast)}
        print("  Sun and Earth from JPL Horizons ...", flush=True)
        rec["sun"] = horizons(10, lat, lon, a.start, stop)
        rec["earth"] = horizons(399, lat, lon, a.start, stop)
        n = min(len(rec["sun"]["el"]), len(rec["earth"]["el"]), a.days * 24)
        for k in ("sun", "earth"):
            rec[k] = {"az": rec[k]["az"][:n], "el": rec[k]["el"][:n]}
        json.dump(rec, open(f"data/{sid}_{a.start}.json", "w"))
        s, e, b = summarise(rec)
        print(f"  sunlit {s:.1f}%   Earth in view {e:.1f}%   both {b:.1f}%   (max horizon {max(rec['horizon']):.2f} deg)")
