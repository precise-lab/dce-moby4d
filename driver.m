close all;
clear all;

addpath("./src");

root_prefix = '/workspace/shared_data/DCE-MOBY4D/vit_z1_r3';
prop_prefix  = 'properties';
moby_anatomy_prefix = fullfile(root_prefix, '../anatomical_structure_body_vessel_flag_0');
moby_lesion_prefix  = fullfile(root_prefix, '../anatomical_structure_lesion_z1_r3');
output_prefix = fullfile(root_prefix, 'phantom');

if(~exist(root_prefix, 'dir'))
    disp(['mkdir ', root_prefix]);
    mkdir(root_prefix);
end

if(~exist(output_prefix, 'dir'))
    disp(['mkdir ', output_prefix]);
    mkdir(output_prefix);
end

% Volume fraction of ICG in blood
f_ICG_in_blood = 50e-6/2e-3;

info = ini2struct(fullfile(moby_anatomy_prefix,'moby.par'));

% Dimension size [voxel]
Nx = str2num(info.array_size);
Ny = str2num(info.array_size);
Nz = str2num(info.endslice)-str2num(info.startslice)+1;
dx = str2num(info.pixel_width)*10; % voxel size in x-dir (mm)
dy = str2num(info.pixel_width)*10; % voxel size in y-dir (mm)
dz = str2num(info.slice_width)*10; % voxel size in z-dir (mm)

% Anatomy frames
anatomy_frame_n = str2num(info.out_frames);
% Time frame
dt      = str2num(info.time_per_frame); % s

fprintf('[Nx, Ny, Nz]: [%d, %d, %d]\n', Nx, Ny, Nz);
fprintf('[dx, dy, dz]: [%.3f, %.3f, %.3f] (mm):\n', dx, dy, dz);
fprintf('Number of anatomical frames: %d\n', anatomy_frame_n);
fprintf('Sampling rate: %.3f s.\n', dt);

% Number of frames
startTime = -2*36; % 2 rotations and 36 seconds per rotation before injection starts
endTime = 17 * 36; %17 rotations and 36 seconds per rotation during and after injection
frame_n = (endTime-startTime)/dt;

% Wavelength
lamb = 800; % [nm]

[ca, c_perf, t] = contrast_agent_curve(startTime, endTime, dt);


% plot and save the CA curve
%figure; plot(t, ca, 'r', 'LineWidth', 3); hold on;
%plot(t, c_perf, 'b--', 'LineWidth', 3); hold on;
%ylim([0, 2.2]);
%ylabel('Concentration (mM)');
%xlabel('Time (s)');
%legend({'C_a(t)', 'C_{perf}(t)'});
%set(gca, 'FontSize', 24);

save(fullfile(output_prefix, 'contrast_agent_curve.mat'), 'ca', 'c_perf', 't' );

properties = load_properties(prop_prefix, f_ICG_in_blood);

fid = fopen(fullfile(root_prefix,'config.ini'),'w');
fprintf(fid, '[path]\n');
fprintf(fid, 'root_folder = %s\n', root_prefix);
fprintf(fid, 'fname_template = moby{0:06d}\n');
fprintf(fid, '[grid]\n');
fprintf(fid, 'frame_n = %d\n', frame_n);
fprintf(fid, 'anatomy_frame_n = %d\n', anatomy_frame_n);
fprintf(fid, 'voxel_size_x = %f\n', dx);
fprintf(fid, 'voxel_size_y = %f\n', dy);
fprintf(fid, 'voxel_size_z = %f\n', dz);
fprintf(fid, 'frame_rate = %f\n', dt);
fprintf(fid, '[labels]\n');
fprintf(fid, 'spleen = %d\n', 37);
fprintf(fid, 'lesion = %d\n', 9);
fprintf(fid, 'intestin = [%d, %d, %d, %d]\n', 38, 39, 41, 42);
fprintf(fid, 'brain = [%d, %d]\n', 43, 81);
fclose(fid);

start_f = 1;
fid = fopen(fullfile(root_prefix,'mu_a_trace'),'w');

for i_frame=start_f:frame_n

    i_anatomy_frame = mod(i_frame-1, anatomy_frame_n)+1;
    label_map = get_label_map(moby_anatomy_prefix, moby_lesion_prefix, [Nx, Ny, Nz], i_anatomy_frame);


    mu_a = assign_mu_a(properties, label_map, lamb, ca(i_frame), c_perf(i_frame));
    mu_sp = assign_mu_sp(properties, label_map, lamb);

    fname = fullfile(output_prefix, ['moby', num2str(i_frame, '%06d'), '.mat']);
    
    label_map = int16(label_map);
    mu_a = single(mu_a);
    mu_sp = single(mu_sp);

    save(fname, 'label_map', 'mu_a', 'mu_sp');

    lesion_volume = sum(label_map==9, 'all')*dx*dy*dz;
    lesion_mu_a = sum(mu_a(label_map==9), 'all')/sum(label_map==9, 'all');
    spleen_volume = sum(label_map==37, 'all')*dx*dy*dz;
    spleen_mu_a = sum(mu_a(label_map==37), 'all')/sum(label_map==37, 'all');
    fprintf(fid, '%06d %1.5e %1.5e %1.5e %1.5e\n', i_frame, lesion_mu_a, ...
                                spleen_mu_a, lesion_volume, spleen_volume);

    fprintf('%06d %1.5e %1.5e %1.5e %1.5e\n', i_frame, lesion_mu_a, ...
                                spleen_mu_a, lesion_volume, spleen_volume);

end

