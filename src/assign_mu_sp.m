% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     May 3, 2022
%

% Dimension size [voxel]
Nx = 372;
Ny = 372;
Nz = 1088;

% Number of frames
frame_n = 10; % 10 frames (respiraotry cycle of 5 sec) repeat

% Wavelength
lamb = 800; % [nm]

% Load label
load(fullfile('properties', 'label.mat'));

% Load optical properties
load(fullfile('properties', 'opt_prop.mat'));

% Reduced scattering coefficient mu_sp [1/mm]
mu_sp_table = mu_s_ref.*(lamb./wavelength_ref).^(-0.7).*(1 - g);

% Assign mu_sp to MOBY phantom with spherical lesion inserted
for frame_i = 1:frame_n
    fprintf('Assigning mu_sp of Frame %d...\n', frame_i);
    
    % Load MOBY phantom with a spherical lesion inserted
    fname_phan = fullfile('anatomical_structure_body_lesion', ...
        ['moby_', num2str(frame_i), '.mat']);
    load(fname_phan, 'phan');
    
    % Assign mu_sp to each tissue
    mu_sp = zeros(Nx, Ny, Nz);
    phan_label = unique(phan);
    for phan_label_i = 1:length(phan_label)
        idx = cellfun(@(x)x == phan_label(phan_label_i), label(:, 3), ...
            'UniformOutput', 1);
        mu_sp(phan == phan_label(phan_label_i)) = mu_sp_table(idx);
    end
    mu_sp = single(mu_sp);
    
    % Save mu_sp map
    fname_mu_sp = fullfile('mu_sp', ...
        ['mu_sp_w', num2str(lamb), '_', num2str(frame_i), '.mat']);
    save(fname_mu_sp, 'mu_sp');
end