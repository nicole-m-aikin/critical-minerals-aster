// ASTER TIR parity recipe — bisbee
// Collection ASTER/AST_L1T_003 (GEE v003). Local pipeline uses LP DAAC v004.
// Clip the study bbox *before* percentiles. Skipping that clip is a known
// local-pipeline bug class; the GEE analogue is reduceRegion without region.
var region = ee.Geometry.Rectangle([-110.15, 31.25, -109.65, 31.7], 'EPSG:4326', false);
var col = ee.ImageCollection('ASTER/AST_L1T_003')
  .filterBounds(region)
  .filterDate('2010-01-01', '2023-12-31')
  .select(['B10', 'B11', 'B12', 'B13', 'B14']);
var image = col
  .map(function(img) {
    var inter = img.geometry().intersection(region, 1);
    return img.set('parity_coverage', inter.area(1).divide(region.area(1)));
  })
  .filter(ee.Filter.gt('parity_coverage', 0.3))
  .sort('parity_coverage', false).first();
image = image.select(['B10', 'B11', 'B12', 'B13', 'B14']).clip(region);
image = image.updateMask(image.select('B10').gt(0));
var silica    = image.select('B13').divide(image.select('B14')).rename('silica');
var carbonate = image.select('B13').divide(image.select('B12')).rename('carbonate');
var mafic     = image.select('B12').divide(image.select('B13')).rename('mafic');
var ratios = silica.addBands(carbonate).addBands(mafic);
var pct = ratios.reduceRegion({
  reducer: ee.Reducer.percentile([70, 90]),
  geometry: region,
  scale: 90.0,
  // crs omitted → GEE defaults to EPSG:4326 (resampling discrepancy)
  maxPixels: 1e9
});
