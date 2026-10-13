"""AVIRIS reflectance download and band-subset I/O.

Supports AVIRIS-5 L2A orthocorrected netCDF (``*_RFL_ORT.nc``) and
AVIRIS-Classic L2 ENVI (``.rfl.hdr`` + ``.rfl.bin``).  Both paths return
``(stack, wavelengths_nm, transform, crs)``.  EarthData credentials come
from ``~/.netrc``, same as the ASTER path.
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
AVIRISC_L2_CONCEPT_ID = "C2711871294-ORNL_CLOUD"
RFL_SUFFIX = "_RFL_ORT.nc"
CLASSIC_BIN_SUFFIXES = (".rfl.bin", "_rfl.bin", ".img")
CLASSIC_HDR_SUFFIXES = (".rfl.hdr", "_rfl.hdr", ".hdr")
NODATA = -9999.0
MAX_WAVELENGTH_DELTA_NM = 8.0

COLLECTION_CONCEPT_IDS = {
    "aviris5_l2a_rfl": AVIRIS5_L2A_CONCEPT_ID,
    "avirisc_l2_rfl": AVIRISC_L2_CONCEPT_ID,
}

COLLECTION_SENSORS = {
    "aviris5_l2a_rfl": "AVIRIS-5",
    "avirisc_l2_rfl": "AVIRIS-Classic",
}

COLLECTION_GSD_M = {
    "aviris5_l2a_rfl": 10.0,
    "avirisc_l2_rfl": 20.0,
}


def wavelength_targets_ok(
    wavelengths_nm,
    targets_nm: tuple[float, ...],
    *,
    max_delta_nm: float = MAX_WAVELENGTH_DELTA_NM,
) -> bool:
    """True when every target has a band within *max_delta_nm*."""
    w = np.asarray(wavelengths_nm, dtype=float)
    if w.size == 0:
        return False
    for t in targets_nm:
        if float(np.min(np.abs(w - float(t)))) > max_delta_nm:
            return False
    return True


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



def _parse_envi_value(raw: str):
    raw = raw.strip()
    if raw.startswith("{") and raw.endswith("}"):
        inner = raw[1:-1].strip()
        if not inner:
            return []
        parts = [p.strip() for p in inner.split(",")]
        out = []
        for p in parts:
            try:
                out.append(float(p))
            except ValueError:
                out.append(p.strip().strip("'\""))
        return out
    try:
        if "." in raw:
            return float(raw)
        return int(raw)
    except ValueError:
        return raw.strip().strip("'\"")


def parse_envi_header(hdr_path: Path) -> dict:
    """Parse an ENVI ``.hdr`` into typed fields for Classic L2 reflectance."""
    hdr_path = Path(hdr_path)
    if not hdr_path.is_file():
        raise FileNotFoundError(f"missing ENVI hdr: {hdr_path}")
    text = hdr_path.read_text(errors="replace")
    if not text.lstrip().upper().startswith("ENVI"):
        raise ValueError(f"not an ENVI header: {hdr_path}")
    # Fold brace-wrapped multi-line values into one line.
    flat = []
    buf = ""
    depth = 0
    for line in text.splitlines():
        if line.strip().upper() == "ENVI":
            continue
        depth += line.count("{") - line.count("}")
        buf = f"{buf} {line.strip()}" if buf else line.strip()
        if depth <= 0 and buf:
            flat.append(buf)
            buf = ""
            depth = 0
    if buf:
        flat.append(buf)
    meta: dict = {}
    for line in flat:
        if "=" not in line:
            continue
        key, val = line.split("=", 1)
        meta[key.strip().lower()] = _parse_envi_value(val)
    required = ("samples", "lines", "bands", "data type", "interleave")
    missing = [k for k in required if k not in meta]
    if missing:
        raise ValueError(f"ENVI header missing {missing}: {hdr_path}")
    waves = meta.get("wavelength") or meta.get("wavelengths") or []
    if not isinstance(waves, list) or not waves:
        raise ValueError(f"ENVI header has no wavelength list: {hdr_path}")
    map_info = meta.get("map info") or []
    transform, crs = _envi_geo_from_map_info(map_info)
    return {
        "samples": int(meta["samples"]),
        "lines": int(meta["lines"]),
        "bands": int(meta["bands"]),
        "data_type": int(meta["data type"]),
        "interleave": str(meta["interleave"]).lower(),
        "header_offset": int(meta.get("header offset") or 0),
        "byte_order": int(meta.get("byte order") or 0),
        "nodata": float(meta.get("data ignore value") or NODATA),
        "wavelengths_nm": [float(x) for x in waves],
        "map_info": map_info,
        "transform": transform,
        "crs": crs,
        "hdr_path": hdr_path,
    }


def _envi_geo_from_map_info(map_info) -> tuple:
    """Pixel-corner affine + CRS from ENVI map info list."""
    from rasterio.transform import Affine

    if not map_info or len(map_info) < 7:
        raise ValueError(f"ENVI map info too short: {map_info!r}")
    # map info = {UTM, 1, 1, easting, northing, xsize, ysize, zone, North, ...}
    ref_x = float(map_info[1])
    ref_y = float(map_info[2])
    easting = float(map_info[3])
    northing = float(map_info[4])
    xsize = float(map_info[5])
    ysize = float(map_info[6])
    # ENVI reference is the center of (ref_x, ref_y) pixel (1-based).
    west = easting - (ref_x - 1.0) * xsize - xsize / 2.0
    north = northing + (ref_y - 1.0) * ysize + ysize / 2.0
    transform = Affine(xsize, 0.0, west, 0.0, -abs(ysize), north)
    crs = "EPSG:4326"
    kind = str(map_info[0]).upper()
    if kind == "UTM" and len(map_info) >= 8:
        zone = int(float(map_info[7]))
        hemi = "N"
        if len(map_info) >= 9:
            hemi = str(map_info[8])[:1].upper()
        epsg = 32600 + zone if hemi != "S" else 32700 + zone
        crs = f"EPSG:{epsg}"
    return transform, crs


def classic_hdr_for_bin(bin_path: Path, hdr_path: Path | None = None) -> Path:
    """Resolve the ENVI header path beside a Classic reflectance binary."""
    bin_path = Path(bin_path)
    if hdr_path is not None:
        hdr_path = Path(hdr_path)
        if not hdr_path.is_file():
            raise FileNotFoundError(f"missing ENVI hdr: {hdr_path}")
        return hdr_path
    candidates = [
        bin_path.with_suffix(".hdr"),
        Path(str(bin_path) + ".hdr"),
        bin_path.with_name(bin_path.name.replace(".bin", ".hdr")),
        bin_path.with_name(bin_path.name.replace(".rfl.bin", ".rfl.hdr")),
    ]
    for cand in candidates:
        if cand.is_file():
            return cand
    raise FileNotFoundError(f"no ENVI hdr beside {bin_path}")


def read_envi_rfl_subset(
    bin_path: Path,
    wavelengths_nm: tuple[float, ...],
    *,
    hdr_path: Path | None = None,
) -> tuple[np.ndarray, np.ndarray, "rasterio.Affine", str]:
    """Read nearest bands from an AVIRIS-Classic ENVI reflectance cube.

    Returns ``(stack, actual_wavelengths, transform, crs)``.
    """
    bin_path = Path(bin_path)
    hdr = parse_envi_header(classic_hdr_for_bin(bin_path, hdr_path))
    if hdr["data_type"] != 4:
        raise ValueError(f"unsupported ENVI data type {hdr['data_type']} (want 4=float32)")
    interleave = str(hdr["interleave"]).lower()
    if interleave not in {"bsq", "bil"}:
        raise ValueError(
            f"unsupported ENVI interleave {hdr['interleave']!r} (want bsq or bil)"
        )
    waves = np.asarray(hdr["wavelengths_nm"], dtype=float)
    if not wavelength_targets_ok(waves, wavelengths_nm):
        raise ValueError(
            f"Classic wavelengths too coarse for targets {wavelengths_nm}; "
            f"nearest deltas exceed {MAX_WAVELENGTH_DELTA_NM} nm"
        )
    indices = [int(np.argmin(np.abs(waves - t))) for t in wavelengths_nm]
    actual = waves[indices]
    lines, samples, bands = hdr["lines"], hdr["samples"], hdr["bands"]
    dtype = "<f4" if hdr["byte_order"] == 0 else ">f4"
    stack = np.empty((len(indices), lines, samples), dtype=np.float32)
    with open(bin_path, "rb") as fh:
        if interleave == "bsq":
            band_bytes = lines * samples * 4
            for out_i, band_i in enumerate(indices):
                fh.seek(hdr["header_offset"] + band_i * band_bytes)
                raw = fh.read(band_bytes)
                if len(raw) != band_bytes:
                    raise ValueError(f"truncated Classic band {band_i} in {bin_path}")
                stack[out_i] = np.frombuffer(raw, dtype=dtype).reshape(lines, samples)
        else:
            # BIL: for each line, bands are contiguous blocks of *samples*.
            line_bytes = bands * samples * 4
            band_row_bytes = samples * 4
            for row in range(lines):
                line_off = hdr["header_offset"] + row * line_bytes
                for out_i, band_i in enumerate(indices):
                    fh.seek(line_off + band_i * band_row_bytes)
                    raw = fh.read(band_row_bytes)
                    if len(raw) != band_row_bytes:
                        raise ValueError(
                            f"truncated Classic BIL band {band_i} row {row} in {bin_path}"
                        )
                    stack[out_i, row, :] = np.frombuffer(raw, dtype=dtype)
    nodata = float(hdr["nodata"])
    stack = np.where(stack == nodata, np.nan, stack)
    return stack, actual, hdr["transform"], hdr["crs"]


def classic_download_urls(granule_ur: str) -> tuple[str, str]:
    """Return ``(hdr_url, bin_url)`` GET DATA links for a Classic granule.

    Classic ORNL granules have no ProducerGranuleId; query CMR by
    ``granule_ur`` (same string stored in the coverage cache).
    """
    query = urllib.parse.urlencode({
        "collection_concept_id": AVIRISC_L2_CONCEPT_ID,
        "granule_ur": granule_ur,
        "page_size": 1,
    })
    payload = _get_json(f"{CMR_GRANULE_URL}?{query}")
    items = payload.get("items") or []
    if not items:
        # Fallback: readable name (some CMR indexes omit exact GranuleUR match).
        query = urllib.parse.urlencode({
            "collection_concept_id": AVIRISC_L2_CONCEPT_ID,
            "readable_granule_name": granule_ur,
            "page_size": 1,
        })
        payload = _get_json(f"{CMR_GRANULE_URL}?{query}")
        items = payload.get("items") or []
    if not items:
        raise FileNotFoundError(f"CMR has no Classic granule {granule_ur!r}")
    hdr_url = bin_url = None
    for rel in items[0].get("umm", {}).get("RelatedUrls") or []:
        if rel.get("Type") != "GET DATA":
            continue
        url = rel.get("URL") or ""
        low = url.lower()
        if low.endswith(".hdr") or low.endswith("_rfl.hdr") or low.endswith(".rfl.hdr"):
            hdr_url = url
        if (
            low.endswith(".bin")
            or low.endswith("_rfl.bin")
            or low.endswith(".rfl.bin")
            or low.endswith(".img")
        ):
            bin_url = url
    if not hdr_url or not bin_url:
        raise FileNotFoundError(
            f"no Classic hdr/bin GET DATA links for {granule_ur} "
            f"(hdr={bool(hdr_url)}, bin={bool(bin_url)})"
        )
    return hdr_url, bin_url


#: Skip Classic cubes larger than this (full-scene orthorectified bins can
#: exceed 100 GB; the SWIR subset only needs a few bands). Override via env
#: ``AVIRIS_CLASSIC_MAX_CUBE_GB``.
CLASSIC_MAX_CUBE_BYTES = int(
    float(__import__("os").environ.get("AVIRIS_CLASSIC_MAX_CUBE_GB", "45")) * 1e9
)


class ClassicCubeTooLarge(RuntimeError):
    """Raised when a Classic reflectance cube exceeds the download size gate."""


def download_classic_rfl(
    granule_ur: str,
    dest_dir: Path,
    session=None,
    *,
    max_cube_bytes: int | None = None,
) -> tuple[Path, Path]:
    """Download Classic ENVI hdr+bin into *dest_dir*; resume if both present."""
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    max_bytes = CLASSIC_MAX_CUBE_BYTES if max_cube_bytes is None else max_cube_bytes
    hdr_url, bin_url = classic_download_urls(granule_ur)
    hdr_dest = dest_dir / Path(urllib.parse.urlparse(hdr_url).path).name
    bin_dest = dest_dir / Path(urllib.parse.urlparse(bin_url).path).name
    if (
        hdr_dest.is_file()
        and bin_dest.is_file()
        and bin_dest.stat().st_size > 10 * 1024 * 1024
    ):
        if bin_dest.stat().st_size > max_bytes:
            raise ClassicCubeTooLarge(
                f"{granule_ur} on-disk bin {bin_dest.stat().st_size/1e9:.1f} GB "
                f"> max {max_bytes/1e9:.1f} GB"
            )
        return hdr_dest, bin_dest

    if session is None:
        import earthaccess

        earthaccess.login(strategy="netrc")
        session = earthaccess.get_requests_https_session()

    for url, dest in ((hdr_url, hdr_dest), (bin_url, bin_dest)):
        if dest.is_file() and (
            dest.suffix.lower() == ".hdr" or dest.stat().st_size > 10 * 1024 * 1024
        ):
            continue
        print(f"  [aviris-classic] downloading {dest.name} …", flush=True)
        with session.get(url, stream=True, timeout=600) as resp:
            resp.raise_for_status()
            cl = resp.headers.get("Content-Length")
            if cl and dest.suffix.lower() == ".bin" and int(cl) > max_bytes:
                raise ClassicCubeTooLarge(
                    f"{granule_ur} Content-Length {int(cl)/1e9:.1f} GB "
                    f"> max {max_bytes/1e9:.1f} GB"
                )
            tmp = dest.with_suffix(dest.suffix + ".part")
            written = 0
            with open(tmp, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=8 * 1024 * 1024):
                    if not chunk:
                        continue
                    written += len(chunk)
                    if dest.suffix.lower() == ".bin" and written > max_bytes:
                        fh.close()
                        tmp.unlink(missing_ok=True)
                        raise ClassicCubeTooLarge(
                            f"{granule_ur} download exceeded max "
                            f"{max_bytes/1e9:.1f} GB (got {written/1e9:.1f} GB)"
                        )
                    fh.write(chunk)
            tmp.replace(dest)
        print(
            f"  [aviris-classic] wrote {dest} "
            f"({dest.stat().st_size / 1e6:.1f} MB)",
            flush=True,
        )
    return hdr_dest, bin_dest


def download_swir_product(
    granule_ur: str,
    dest_dir: Path,
    *,
    collection_key: str,
    session=None,
):
    """Download AVIRIS-5 netCDF or Classic ENVI pair for *collection_key*."""
    if collection_key == "aviris5_l2a_rfl":
        return download_rfl_ort(granule_ur, dest_dir, session=session)
    if collection_key == "avirisc_l2_rfl":
        return download_classic_rfl(granule_ur, dest_dir, session=session)
    raise ValueError(f"unsupported SWIR collection {collection_key!r}")


def read_swir_subset(
    path: Path,
    wavelengths_nm: tuple[float, ...],
    *,
    collection_key: str,
    hdr_path: Path | None = None,
):
    """Collection-aware reflectance band subset."""
    if collection_key == "aviris5_l2a_rfl":
        stack, actual, transform, crs = read_rfl_subset(path, wavelengths_nm)
        if not wavelength_targets_ok(actual, wavelengths_nm):
            raise ValueError(
                f"AVIRIS-5 wavelengths too far from targets {wavelengths_nm}"
            )
        return stack, actual, transform, crs
    if collection_key == "avirisc_l2_rfl":
        return read_envi_rfl_subset(path, wavelengths_nm, hdr_path=hdr_path)
    raise ValueError(f"unsupported SWIR collection {collection_key!r}")


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
    # Per must-contain point, seed at most ONE granule (largest intersection with
    # the site). Seeding every overlapping Classic strip can force 5–10×20 GB
    # downloads for a single town point.
    points = [Point(lon, lat) for lon, lat in must_contain]
    seed_best: dict[int, tuple[str, object, float]] = {}
    for rec in recs:
        wkt = rec.get("geometry_wkt")
        if not wkt:
            continue
        geom = shapely_wkt.loads(wkt)
        inter = geom.intersection(site)
        if inter.is_empty:
            continue
        ur = rec["granule_ur"]
        inter_eq = shp_transform(to_eq, inter)
        candidates.append((ur, inter_eq))
        for pi, pt in enumerate(points):
            if not geom.contains(pt):
                continue
            area = float(inter_eq.area)
            prev = seed_best.get(pi)
            if prev is None or area > prev[2]:
                seed_best[pi] = (ur, inter_eq, area)

    selected: list[str] = []
    union = None
    remaining = list(candidates)
    # Deduplicate seeds that share a granule across points
    for ur, geom, _area in seed_best.values():
        if ur in selected:
            continue
        match = next((item for item in remaining if item[0] == ur), None)
        if match is None:
            continue
        remaining = [item for item in remaining if item[0] != ur]
        selected.append(ur)
        union = geom if union is None else union.union(geom)

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
