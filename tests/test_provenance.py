import hashlib

from critical_minerals_aster.config import ClassificationParams, SiteConfig
from critical_minerals_aster.paths import SitePaths
from critical_minerals_aster.pipeline import input_provenance

MRDS_CSV = (
    "dep_id,site_name,latitude,longitude,dev_stat,score,commod1\n"
    "1,Pit,40.0,-115.0,Producer,B,Gold\n"
    "2,Mill,40.1,-115.1,Plant,C,Gold\n"
    "3,Pit,40.0,-115.0,Occurrence,D,Gold\n"
)


def _site(tmp_path):
    (tmp_path / "sites").mkdir()
    (tmp_path / "sites" / "demo.yaml").write_text("id: demo\n")
    site = SiteConfig(
        id="demo",
        name="Demo",
        bbox_wgs84=(-116.0, 39.0, -114.0, 41.0),
        layout="nested",
        classification=ClassificationParams(low_pct=65.0, high_pct=85.0, strong_score_min=2),
    )
    return site, SitePaths(repo_root=tmp_path, site=site)


def test_input_provenance_records_hashes_params_and_validation(tmp_path):
    site, paths = _site(tmp_path)
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "mrds.csv").write_text(MRDS_CSV)

    inputs = input_provenance(site, paths, "AST_L1T_00401012010052826_20250701203557")

    assert inputs["site_config_sha256"] == hashlib.sha256(b"id: demo\n").hexdigest()
    assert inputs["mrds_sha256"] == hashlib.sha256(MRDS_CSV.encode()).hexdigest()
    assert inputs["classification"] == {"low_pct": 65.0, "high_pct": 85.0, "strong_score_min": 2}
    assert inputs["aster_product"] == "AST_L1T"
    assert inputs["aster_product_version"] == "004"
    report = inputs["mrds_validation_national"]
    assert report["n_input"] == 3
    assert report["non_deposit_plant"] == 1
    assert report["duplicate_record"] == 1
    assert report["n_output"] == 1


def test_input_provenance_without_mrds_or_granule(tmp_path):
    site, paths = _site(tmp_path)
    inputs = input_provenance(site, paths, "demo_mosaic")
    assert inputs["mrds_sha256"] is None
    assert inputs["mrds_validation_national"] is None
    assert inputs["aster_product_version"] is None
    assert inputs["structure_source"] == []
