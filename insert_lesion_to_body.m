% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     May 3, 2022
%

% Dimension size [voxel]
Nx = 372;
Ny = 372;
Nz = 1088;

% Number of frames
frame_n = 10;

% Insert a spherical lesion to MOBY body
for frame_i = 1:frame_n
    fprintf('Inserting a lesion to body of Frame %d...\n', frame_i);
    
    % Load MOBY phantom and spherical lesion phantom
    fname_body = fullfile('anatomical_structure_body', ...
        ['moby_act_', num2str(frame_i), '.bin']);
    fname_lesn = fullfile('anatomical_structure_lesion', ...
        ['lesion_act_', num2str(frame_i), '.bin']);
    
    fid = fopen(fname_body, 'rb'); phan = fread(fid, 'float'); fclose(fid);
    fid = fopen(fname_lesn, 'rb'); lesn = fread(fid, 'float'); fclose(fid);
    
    % Insert the spherical lesion into the phantom
    phan(lesn ~= 0) = lesn(lesn ~= 0);
    phan = reshape(phan, [Nx, Ny, Nz]);
    
    % Save mu_sp map
    fname_phan = fullfile('anatomical_structure_body_lesion', ...
        ['moby_', num2str(frame_i), '.mat']);
    save(fname_phan, 'phan');
end