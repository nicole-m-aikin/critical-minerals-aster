"""Dry-run granule inventory for mountain_pass — checks size/coverage filters."""
import earthaccess
earthaccess.login(strategy="netrc")

from critical_minerals_aster.spectral import score_granule

bbox = (-115.75, 35.28, -115.30, 35.68)
results = earthaccess.search_data(
    short_name="AST_L1T",
    bounding_box=bbox,
    temporal=("2010-01-01", "2023-12-31"),
    count=20,
)
print(f"Found {len(results)} granules\n")
print("%-8s %-6s %-6s  %s" % ("MB", "cov", "bands", "granule"))
print("-" * 75)
for g in results:
    try:
        cov, _, bands = score_granule(g, bbox)
        sz = g.size()
        gid = g["umm"]["GranuleUR"]
        flags = []
        if sz >= 20:
            flags.append("LARGE")
        if cov <= 0.30:
            flags.append("LOW-COV")
        note = "  <- " + ", ".join(flags) if flags else ""
        print("%-8.1f %-6.2f %-6d  %s%s" % (sz, cov, bands, gid, note))
    except Exception as e:
        print(f"  (error: {e})")
