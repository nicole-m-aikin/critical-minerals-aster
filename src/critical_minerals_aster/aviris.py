"""AVIRIS-5 L2A reflectance download and band-subset I/O.

Used by the Goldfield SWIR pilot.  Downloads only the orthocorrected
reflectance netCDF (``*_RFL_ORT.nc``, ~2 GB); the matching uncertainty
cube is not fetched.  EarthData credentials come from ``~/.netrc``, same
as the ASTER path.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import Affine, from_origin

from critical_minerals_aster.hyperspectral import (
    CMR_GRANULE_URL,
    _get_json,
)

AVIRIS5_L2A_CONCEPT_ID = "C4184090548-ORNL_CLOUD"
RFL_SUFFIX = "_RFL_ORT.nc"
NODATA = -9999.0


def rfl_download_url(granule_ur: str) -> str:
    """HTTPS URL of the ``*_RFL_ORT.nc`` file for *granule_ur*."""
    query = urllib.parse.urlencode({
        "collection_concept_id": AVIRIS5_L2A_CONCEPT_ID,
        "producer_granule_id": granule_ur,
        "page_size": 1,
    })
    payload = _get_json(f"{CMR_GRANULE_URL}?{query}")
    items = payload.get("items") or []
    if not items:
        raise FileNotFoundError(f"CMR has no granule {granule_ur!r}")
    for rel in items[0].get("umm", {}).get("RelatedUrls") or []:
        url = rel.get("URL") or ""
        if rel.get("Type") == "GET DATA" and url.endswith(RFL_SUFFIX):
            return url
    raise FileNotFoundError(f"no {RFL_SUFFIX} GET DATA link for {granule_ur}")


def download_rfl_ort(
    granule_ur: str,
    dest_dir: Path,
    session=None,
) -> Path:
    """Download the reflectance netCDF for *granule_ur* into *dest_dir*.

    Skips the transfer when the file is already present and larger than
    100 MB (a complete AVIRIS-5 scene is ~2 GB; anything smaller is a
    truncated leftover).
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    url = rfl_download_url(granule_ur)
    dest = dest_dir / Path(urllib.parse.urlparse(url).path).name
    if dest.is_file() and dest.stat().st_size > 100 * 1024 * 1024:
        return dest

    if session is None:
        import earthaccess

        earthaccess.login(strategy="netrc")
        session = earthaccess.get_requests_https_session()

    print(f"  [aviris] downloading {dest.name} …", flush=True)
    with session.get(url, stream=True, timeout=600) as resp:
        resp.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=8 * 1024 * 1024):
                if chunk:
                    fh.write(chunk)
        tmp.replace(dest)
    print(f"  [aviris] wrote {dest} ({dest.stat().st_size / 1e9:.2f} GB)", flush=True)
    return dest


def read_rfl_subset(
    nc_path: Path,
    wavelengths_nm: tuple[float, ...],
) -> tuple[np.ndarray, np.ndarray, rasterio.Affine, str]:
    """Read the bands nearest *wavelengths_nm* from an AVIRIS-5 RFL_ORT file.

    Returns ``(stack, actual_wavelengths, transform, crs)`` where *stack*
    is ``(band, y, x)`` float32 with NaNs at nodata.
    """
    import netCDF4

    ds = netCDF4.Dataset(str(nc_path), "r")
    try:
        # ORNL's RFL_ORT files nest wavelength + reflectance under a
        # ``reflectance`` group (the user-guide table shows them at root).
        grp = ds.groups["reflectance"] if "reflectance" in ds.groups else ds
        w = np.asarray(grp["wavelength"][:], dtype=float)
        rfl = grp["reflectance"]
        band_axis = int(np.argmin(np.abs(np.array(rfl.shape) - w.size)))
        if rfl.shape[band_axis] != w.size:
            raise ValueError(
                f"cannot match wavelength axis in shape {rfl.shape} vs {w.size}"
            )

        indices = [int(np.argmin(np.abs(w - t))) for t in wavelengths_nm]
        actual = w[indices]
        # Read one band at a time so we never materialise the full cube.
        slices = []
        for i in indices:
            indexer = [slice(None)] * rfl.ndim
            indexer[band_axis] = i
            band = np.asarray(rfl[tuple(indexer)], dtype=np.float32)
            if band.ndim != 2:
                raise ValueError(f"expected 2-D band, got {band.shape}")
            slices.append(band)
        stack = np.stack(slices, axis=0)
        stack = np.where(stack == NODATA, np.nan, stack)

        easting = np.asarray(ds["easting"][:], dtype=float)
        northing = np.asarray(ds["northing"][:], dtype=float)
        # from_origin expects a north-up raster (row 0 = north). Flip if the
        # netCDF axis runs south → north.
        if northing.size >= 2 and (northing[1] - northing[0]) > 0:
            stack = stack[:, ::-1, :]
            northing = northing[::-1]
        transform, crs = _geo_from_axes(easting, northing, ds)
    finally:
        ds.close()
    return stack, actual, transform, crs


def _geo_from_axes(easting, northing, ds) -> tuple[rasterio.Affine, str]:
    """Pixel-corner affine and CRS from the 1-D UTM axes."""
    if easting.size < 2 or northing.size < 2:
        raise ValueError("easting/northing axes are too short to georeference")
    dx = float(easting[1] - easting[0])
    dy = float(northing[1] - northing[0])
    west = float(easting[0]) - dx / 2.0
    north = float(northing[0]) - dy / 2.0
    transform = from_origin(west, north, abs(dx), abs(dy))

    tm = ds["transverse_mercator"]
    crs = None
    for key in ("spatial_ref", "crs_wkt"):
        if hasattr(tm, key):
            crs = str(getattr(tm, key))
            break
    if crs is None:
        lon0 = float(getattr(tm, "longitude_of_central_meridian", -117.0))
        zone = int((lon0 + 183) / 6)
        crs = f"EPSG:{32600 + zone}"
    return transform, crs


def write_geotiff(
    path: Path,
    array: np.ndarray,
    transform: rasterio.Affine,
    crs: str,
    nodata: float = np.nan,
) -> Path:
    """Write a single-band float32 GeoTIFF."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = np.asarray(array, dtype=np.float32)
    if data.ndim != 2:
        raise ValueError(f"expected 2-D array, got {data.shape}")
    with rasterio.open(
        path, "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype="float32",
        crs=crs,
        transform=transform,
        nodata=-9999.0,
        compress="deflate",
    ) as dst:
        out = np.where(np.isfinite(data), data, -9999.0)
        dst.write(out, 1)
    return path


def list_granules_covering(
    cache_path: Path,
    lon: float,
    lat: float,
) -> list[str]:
    """Cached granule_ur values whose footprint contains (lon, lat)."""
    from shapely.geometry import Point
    from shapely import wkt as shapely_wkt

    recs = json.loads(Path(cache_path).read_text())
    pt = Point(lon, lat)
    hits = []
    for rec in recs:
        wkt = rec.get("geometry_wkt")
        if not wkt:
            continue
        if shapely_wkt.loads(wkt).contains(pt):
            hits.append(rec["granule_ur"])
    return sorted(set(hits))


def select_granules_painting_bbox(
    cache_path: Path,
    bbox_wgs84: tuple[float, float, float, float],
    must_contain: tuple[tuple[float, float], ...] = (),
    min_gain: float = 0.02,
) -> list[str]:
    """Pick a small set of flight lines that together paint *bbox_wgs84*.

    ASTER's mosaic gate is ``coverage > 30%`` per granule.  An AVIRIS-5 swath
    is a thin strip, so no single line clears that bar.  The analogous rule
    is a greedy set-cover: keep adding the footprint that covers the most
    *new* bbox area until the next line would add less than *min_gain*.

    *must_contain* points (town, district) are seeded first so the science
    target is never dropped for a slightly larger empty strip.
    """
    import pyproj
    from shapely.geometry import Point, box
    from shapely import wkt as shapely_wkt
    from shapely.ops import transform as shp_transform

    recs = json.loads(Path(cache_path).read_text())
    site = box(*bbox_wgs84)
    to_eq = pyproj.Transformer.from_crs(
        "EPSG:4326", "EPSG:5070", always_xy=True,
    ).transform
    site_eq = shp_transform(to_eq, site)
    site_area = site_eq.area
    if site_area <= 0:
        return []

    candidates: list[tuple[str, object]] = []
    seeds: list[str] = []
    points = [Point(lon, lat) for lon, lat in must_contain]
    for rec in recs:
        wkt = rec.get("geometry_wkt")
        if not wkt:
            continue
        geom = shapely_wkt.loads(wkt)
        inter = geom.intersection(site)
        if inter.is_empty:
            continue
        ur = rec["granule_ur"]
        candidates.append((ur, shp_transform(to_eq, inter)))
        if any(geom.contains(pt) for pt in points):
            seeds.append(ur)

    selected: list[str] = []
    union = None
    remaining = list(candidates)
    for ur in seeds:
        match = next((item for item in remaining if item[0] == ur), None)
        if match is None:
            continue
        remaining = [item for item in remaining if item[0] != ur]
        selected.append(ur)
        union = match[1] if union is None else union.union(match[1])

    while remaining:
        best_i = None
        best_gain = 0.0
        best_union = None
        base = 0.0 if union is None else union.area
        for i, (_ur, geom) in enumerate(remaining):
            nxt = geom if union is None else union.union(geom)
            gain = (nxt.area - base) / site_area
            if gain > best_gain:
                best_i, best_gain, best_union = i, gain, nxt
        if best_i is None or best_gain < min_gain:
            break
        selected.append(remaining.pop(best_i)[0])
        union = best_union
    return selected


def pick_granule_covering(
    cache_path: Path,
    lon: float,
    lat: float,
) -> str:
    """Return the first cached granule_ur whose footprint contains (lon, lat)."""
    hits = list_granules_covering(cache_path, lon, lat)
    if not hits:
        raise FileNotFoundError(
            f"no cached footprint contains ({lon}, {lat}) in {cache_path}"
        )
    return hits[0]


def mosaic_to_bbox(
    arrays: list[np.ndarray],
    transforms: list[rasterio.Affine],
    crs,
    bbox_wgs84: tuple[float, float, float, float],
    resolution: float = 10.0,
) -> tuple[np.ndarray, rasterio.Affine, str]:
    """Warp one or more 2-D rasters onto *bbox_wgs84* in the source CRS.

    First finite value wins when scenes overlap.  Returns
    ``(mosaic, transform, crs_string)`` with NaNs outside every swath.
    """
    from rasterio.crs import CRS
    from rasterio.transform import from_origin
    from rasterio.warp import Resampling, reproject, transform_bounds

    if not arrays:
        raise ValueError("mosaic_to_bbox needs at least one array")
    dst_crs = CRS.from_user_input(crs)
    west, south, east, north = transform_bounds(
        "EPSG:4326", dst_crs, *bbox_wgs84, densify_pts=21,
    )
    width = max(1, int(np.ceil((east - west) / resolution)))
    height = max(1, int(np.ceil((north - south) / resolution)))
    dst_transform = from_origin(west, north, resolution, resolution)
    mosaic = np.full((height, width), np.nan, dtype=np.float32)
    for arr, src_t in zip(arrays, transforms):
        tmp = np.full((height, width), np.nan, dtype=np.float32)
        reproject(
            source=np.asarray(arr, dtype=np.float32),
            destination=tmp,
            src_transform=src_t,
            src_crs=dst_crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            src_nodata=np.nan,
            dst_nodata=np.nan,
            resampling=Resampling.nearest,
        )
        mosaic = np.where(np.isfinite(mosaic), mosaic, tmp)
    return mosaic, dst_transform, dst_crs.to_string()


def valid_window(array: np.ndarray, pad: int = 8) -> tuple[int, int, int, int]:
    """Return (row0, row1, col0, col1) covering finite pixels, plus *pad*."""
    work = array[0] if array.ndim == 3 else array
    valid = np.isfinite(work)
    if not valid.any():
        return 0, work.shape[0], 0, work.shape[1]
    rows = np.where(valid.any(axis=1))[0]
    cols = np.where(valid.any(axis=0))[0]
    r0 = max(0, int(rows[0]) - pad)
    r1 = min(work.shape[0], int(rows[-1]) + 1 + pad)
    c0 = max(0, int(cols[0]) - pad)
    c1 = min(work.shape[1], int(cols[-1]) + 1 + pad)
    return r0, r1, c0, c1


def crop_to_valid(array: np.ndarray, transform: rasterio.Affine, pad: int = 8):
    """Trim empty (all-NaN) margins so the flight-line fills the frame."""
    r0, r1, c0, c1 = valid_window(array, pad=pad)
    new_t = transform * Affine.translation(c0, r0)
    if array.ndim == 3:
        return array[:, r0:r1, c0:c1], new_t
    return array[r0:r1, c0:c1], new_t
