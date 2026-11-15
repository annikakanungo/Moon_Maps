# Lunar Landing Site Navigator

NASA Space Apps Challenge 2026 project for the challenge **Lunar Landing Site Navigator**.

A tool for comparing lunar south pole landing sites and dates by where the Sun and Earth sit above the local horizon, which tells you about power potential and direct-to-Earth communication windows.

**Live demo:** https://annikakanungo.github.io/Moon_Maps/

## What it does
- Compare two scenarios side by side (site + start date, 14 / 30 / 60 day window).
- Charts of Sun and Earth elevation over time, with a strip showing when both are available.
- Metrics: percent sunlit, percent Earth in view, percent both, longest dark spell, longest Earth gap, peak Sun elevation.
- Built-in sites: Shackleton Crater (point near the pole), Malapert Massif, Amundsen Rim (site AR01), Nobile Rim 2 (study site), plus Apollo 11 as an equatorial sanity check. Custom coordinates and a horizon mask are supported.

## How it works
1. **Web app (`index.html`)** computes Sun and Earth positions in the browser from an analytic lunar ephemeris (Meeus, truncated series) plus lunar libration. Checked against a published worked example and against astropy: Moon position within 0.03 degrees and 55 km. The horizon is flat unless the user sets a horizon mask.
2. **Data pipeline (`prep_site_data.py`)** builds terrain-aware data from NASA sources:
   - Terrain horizon for each azimuth, ray-cast over the LRO LOLA 80 m south pole elevation map.
   - Hourly Sun and Earth azimuth and elevation from NASA JPL Horizons.
   - Prints percent sunlit, Earth in view, and both with the real horizon.

### Example run (30 days from 2026-11-14, 2 m mast, 30 km search radius)
| Site | Sunlit | Earth in view | Both |
|---|---|---|---|
| Malapert Massif | 72.6% | 100.0% | 72.6% |
| Amundsen Rim (AR01) | 53.3% | 48.1% | 31.8% |
| Shackleton (point near pole) | 48.6% | 56.2% | 39.4% |
| Nobile Rim 2 | 41.5% | 21.8% | 6.2% |

The data files from this run are in this repository.

## Run the pipeline yourself
```
pip install numpy rasterio pyproj
mkdir data
curl -L -o data/LDEM_80S_80MPP_ADJ.TIF https://pgda.gsfc.nasa.gov/data/LOLA_20mpp/LDEM_80S_80MPP_ADJ.TIF
python3 prep_site_data.py --start 2026-11-14 --days 30
```
The terrain file is about 181 MB and is not stored in this repository.

## Limitations
- Sun counts as up when any part of its disc clears the horizon; Earth counts as up when its centre does. No antenna pointing limits.
- Sites near the pole are very sensitive to the exact coordinates; the Shackleton entry is a point about 3 km from the pole, not a validated crater-floor or rim site.
- The terrain grid is 80 m per pixel; smaller features are not resolved.
- The Amundsen Rim and Nobile Rim 2 coordinates come from published landing-site studies, not official NASA selections.
- Not an official NASA product.

## Data and references
- LOLA south pole elevation products: Barker et al. (2023), Planetary Science Journal 4, 183. NASA Goddard PGDA, data DOI 10.60903/gsfcpgda-lola-spole.
- NASA JPL Horizons system (ssd.jpl.nasa.gov).
- Meeus, J., Astronomical Algorithms (lunar position and libration).
