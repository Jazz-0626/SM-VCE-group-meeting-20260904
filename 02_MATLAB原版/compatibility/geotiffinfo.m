function info = geotiffinfo(filename)
%GEOTIFFINFO Compatibility shim exposing the RefMatrix used by SM-VCE.
raw = imfinfo(filename);
if ~isfield(raw, 'ModelPixelScaleTag') || ~isfield(raw, 'ModelTiepointTag')
    error('GeoTIFF metadata missing ModelPixelScaleTag or ModelTiepointTag: %s', filename);
end
scale = double(raw.ModelPixelScaleTag(:));
tie = double(raw.ModelTiepointTag(:));
if numel(scale) < 2 || numel(tie) < 6
    error('Unsupported GeoTIFF georeferencing tags: %s', filename);
end
% Mapping Toolbox RefMatrix convention used by the original code:
% [0 dLat; dLon 0; lon0 lat0]. North-up rasters have negative dLat.
info = raw;
info.RefMatrix = [0, -abs(scale(2)); abs(scale(1)), 0; tie(4), tie(5)];
end
