function result_file = run_groupmeeting_experiment(observation_mode, roi, output_dir)
%RUN_GROUPMEETING_EXPERIMENT Reproducible SM-VCE run using Liu's original functions.
%
% observation_mode: 'original7' excludes observations 3 and 6; 'all9' uses all.
% roi: [] for full image, or [row_start row_end col_start col_end] (1-based).
% output_dir: directory receiving the MAT result and run metadata.

if nargin < 1 || isempty(observation_mode), observation_mode = 'original7'; end
if nargin < 2, roi = []; end
if nargin < 3 || isempty(output_dir), output_dir = fullfile(pwd, 'results'); end

runner_dir = fileparts(mfilename('fullpath'));
demo_dir = fullfile(runner_dir, 'SMVCE_MATLAB_原版');
compat_dir = fullfile(runner_dir, 'compatibility');
addpath(compat_dir);
addpath(fullfile(demo_dir, 'SMVCE_code'));
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(demo_dir);

[data,inc,azi,losazienu,leftorright,coor,dem,mask,fault,datainfo] = SMVCE_readdispdata();
data(data == 0) = NaN;

switch lower(observation_mode)
    case 'original7'
        excluded = [3, 6];
    case 'all9'
        excluded = [];
    otherwise
        error('observation_mode must be original7 or all9');
end

if ~isempty(excluded)
    data(:,:,excluded) = [];
    inc(:,:,excluded) = [];
    azi(:,:,excluded) = [];
    losazienu(excluded) = [];
    leftorright(excluded) = [];
    datainfo(excluded,:) = [];
end

if ~isempty(roi)
    assert(numel(roi) == 4, 'roi must be [row_start row_end col_start col_end]');
    r1 = roi(1); r2 = roi(2); c1 = roi(3); c2 = roi(4);
    data = data(r1:r2,c1:c2,:);
    inc = inc(r1:r2,c1:c2,:);
    azi = azi(r1:r2,c1:c2,:);
    dem = dem(r1:r2,c1:c2);
    mask = mask(r1:r2,c1:c2);
    coor.corner_lon = coor.corner_lon + (c1 - 1) * coor.post_lon;
    coor.corner_lat = coor.corner_lat + (r1 - 1) * coor.post_lat;
    coor.nlines = size(data,1);
    coor.width = size(data,2);
end

windowsize = 41;
flag_if_2D = 0;
fsmpara = 2;
flag_interWeight = 0;
flag_adpws = 0;
flag_smad = 0;

if ~exist(output_dir, 'dir'), mkdir(output_dir); end
data_file = fullfile(output_dir, sprintf('DATA_%s.mat', observation_mode));
save(data_file, 'data', 'inc', 'azi', 'losazienu', 'flag_if_2D', ...
    'leftorright', 'coor', 'datainfo', 'dem', 'mask', 'fault', ...
    'windowsize', 'fsmpara', 'flag_interWeight', 'flag_adpws', 'flag_smad', '-v7.3');

fprintf('Running MATLAB SM-VCE: mode=%s, size=%dx%d, observations=%d\n', ...
    observation_mode, size(data,1), size(data,2), size(data,3));
tic;
Result_smvce = SMVCE_solve3D(data_file);
elapsed_seconds = toc;

enu = Result_smvce.enu;
enu_std = Result_smvce.var.enu;
sita = Result_smvce.sita;
result_file = fullfile(output_dir, sprintf('MATLAB_%s_result.mat', observation_mode));
save(result_file, 'enu', 'enu_std', 'sita', 'coor', 'datainfo', ...
    'elapsed_seconds', 'observation_mode', 'roi', '-v7');

fid = fopen(fullfile(output_dir, sprintf('MATLAB_%s_run.txt', observation_mode)), 'w');
fprintf(fid, 'implementation=MATLAB original functions\n');
fprintf(fid, 'observation_mode=%s\n', observation_mode);
fprintf(fid, 'excluded_observations=%s\n', mat2str(excluded));
fprintf(fid, 'windowsize=%d\n', windowsize);
fprintf(fid, 'fsmpara=%d\n', fsmpara);
fprintf(fid, 'roi=%s\n', mat2str(roi));
fprintf(fid, 'rows=%d\ncols=%d\nobservations=%d\n', size(data,1), size(data,2), size(data,3));
fprintf(fid, 'elapsed_seconds=%.6f\n', elapsed_seconds);
fclose(fid);
fprintf('Saved MATLAB result: %s\n', result_file);
end
