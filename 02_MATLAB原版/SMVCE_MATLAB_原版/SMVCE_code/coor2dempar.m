function []=coor2dempar(coor,pout)
%
% coor.corner_lon;
% coor.corner_lat;
% coor.post_lat;
% coor.post_lat;
% coor.width;
% coor.nlines;
%coor2dempar
%--By Liu Jihong 2021/09/26 @ Central South University
%--This function is used to write coor into a file

fid=fopen(pout,'w');
fprintf(fid,'Gamma DIFF&GEO DEM/MAP parameter file\n');
fprintf(fid,'title: DEM\n');
fprintf(fid,'DEM_projection:     EQA\n');
fprintf(fid,'data_format:        REAL*4\n');
fprintf(fid,'DEM_hgt_offset:          0.00000\n');
fprintf(fid,'DEM_scale:               1.00000\n');
fprintf(fid,'width:                %d\n',coor.width);
fprintf(fid,'nlines:               %d\n',coor.nlines);
fprintf(fid,'corner_lat:     %.7f  decimal degrees\n',coor.corner_lat);
fprintf(fid,'corner_lon:   %.7f  decimal degrees\n',coor.corner_lon);
fprintf(fid,'post_lat:   %.7f  decimal degrees\n',coor.post_lat);
fprintf(fid,'post_lon:    %.7f  decimal degrees\n',coor.post_lon);
fprintf(fid,'\n');
fprintf(fid,'ellipsoid_name: WGS 84\n');
fprintf(fid,'ellipsoid_ra:        6378137.000   m\n');
fprintf(fid,'ellipsoid_reciprocal_flattening:  298.2572236\n');
fprintf(fid,'\n');
fprintf(fid,'datum_name: WGS 1984\n');
fprintf(fid,'datum_shift_dx:              0.000   m\n');
fprintf(fid,'datum_shift_dy:              0.000   m\n');
fprintf(fid,'datum_shift_dz:              0.000   m\n');
fprintf(fid,'datum_scale_m:         0.00000e+00\n');
fprintf(fid,'datum_rotation_alpha:  0.00000e+00   arc-sec\n');
fprintf(fid,'datum_rotation_beta:   0.00000e+00   arc-sec\n');
fprintf(fid,'datum_rotation_gamma:  0.00000e+00   arc-sec\n');
fprintf(fid,'datum_country_list: Global Definition, WGS84, World\n');
fclose(fid);
end