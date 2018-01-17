% close all;

addpath('./1. BasicModules/kd_tree');
addpath('./1. BasicModules');
addpath('./2. Non_rigid_registration');
%  distcomp.feature( 'LocalUseMpiexec', false )
%% Example
if ~exist('VS','var')
    [VS, FS, NS] = read_obj_file('cup_source_nohandle.obj');
%     [VS2, FS2, NS2] = read_obj_file('face-09-surprise.obj');
    [VT, FT, NT] = read_obj_file('cup_target_nohandle.obj');    
end

load('Cup_FWland.mat');
load('Cup_land.mat');

marker =  [FW_land landmark];
% marker = init_marker(VS, FS, VT, FT, 'Face_Marker.mat');
wl = 1;
wc = 0.5;
dist = 0.1;
theta = 45;
x = laplacian_registration(VS, FS, NS, VT, FT, NT, marker, wl, wc, dist, theta);
% [ VS_Reg, VT_Reg ] = non_rigid_registration(VS, FS, VT, FT, 1.0, 0.01, [1 500 3000 5000], marker, 'DF_reg_phase2.mat');
dlacorres = build_correspondence(VS_Reg, FS, VT_Reg, FT, 10, 0.05, 'Face_ICIP_corres.mat');
[ x, nx ] = deformation_transfer(VS, FS, VT, FT, VS2, FS2, corres);
write_obj_file(x, FT, nx, 'head-09-surprise.obj');

fprintf('End of demo..\n');
system('pause');

clear VS VT S_factor T_factor FS FT NS NT maker;