function data = geotiffread(filename)
%GEOTIFFREAD Compatibility shim for MATLAB installations without Mapping Toolbox.
% The SM-VCE demo only needs the raster values; IMREAD preserves the
% single-precision samples and NaN nodata values in the supplied GeoTIFFs.
data = imread(filename);
end
