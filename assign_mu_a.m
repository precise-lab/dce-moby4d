% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     May 3, 2022
%

% Dimension size [voxel]
Nx = 372;
Ny = 372;
Nz = 1088;

% Anatomy frames
anatomy_frame_number = 277;

% Number of frames
endTime = 17 * 36; %17 rotations and 36 seconds per rotation
dt      = 0.1; % s
frame_n = endTime/dt;

% Wavelength
lamb = 800; % [nm]

% Load label
load(fullfile('properties', 'label.mat'));

% Load constants
load(fullfile('properties', 'constants.mat'), ...
    'wavelength', 'e_hbo2', 'e_hb', 'mu_a_w');

% Load molar extinction coefficient of ICG
load(fullfile('properties', 'constant_icg.mat'), 'e_icg');
e_icg = e_icg(:, 1); % 65 uM

% Load functional properties
load(fullfile('properties', 'func_prop.mat'));

% Load contrast_agent_curve
ca, cperf = contrast_agent_curve(endTime, dt);

% Volume fraction of ICG in blood
f_ICG_in_blood = 50e-6/2e-3;

% Optical absorption coefficient mu_a [1/mm] excluding contrast agent
mu_a_table = f_b*log(10)*c_thb_b ...
    .*(s*e_hbo2(wavelength == lamb) + (1 - s)*e_hb(wavelength == lamb)) ...
    + f_w*mu_a_w(wavelength == lamb);

% Optical absorption coefficient mu_a [1/mm] of contrast agent without
% its concentration
mu_a_CA = log(10)*e_icg(wavelength == lamb);

% Assign mu_a to MOBY phantom with spherical lesion inserted
for frame_i = 1:frame_n
    fprintf('Assigning mu_a of Frame %d...\n', frame_i);
    
    % Load MOBY phantom with a spherical lesion inserted
    frame_load_i = mod(frame_i, anatomy_frame_number);
    if frame_load_i == 0
        fname_phan = fullfile('anatomical_structure_body_lesion', ...
            ['moby_1', num2str(frame_load_i), '.mat']);
    else
        fname_phan = fullfile('anatomical_structure_body_lesion', ...
            ['moby_', num2str(frame_load_i), '.mat']);
    end
    load(fname_phan, 'phan');
    
    % Update concentration of contrast agent C_CA
    C_CA = f_ICG_in_blood*f_b*ca(frame_i)*1e-3;
    C_CA(strcmp(label(:,1), 'lesn')) =C_CA(strcmp(label(:,1), 'lesn')) + ...
        f_ICG_in_blood*(1. - f_w(strcmp(label(:,1), 'lesn'))) ...
        *cperf(frame_i)*1e-3;
    
    % Assign mu_a to each tissue
    mu_a = zeros(Nx, Ny, Nz);
    phan_label = unique(phan);
    for phan_label_i = 1:length(phan_label)
        idx = cellfun(@(x)x == phan_label(phan_label_i), label(:, 3), ...
            'UniformOutput', 1);
        mu_a(phan == phan_label(phan_label_i)) ...
            = mu_a_table(idx) + C_CA(idx)*mu_a_CA;
    end
    mu_a = single(mu_a);
    
    % Save mu_a map
    fname_mu_a = fullfile('mu_a', ...
        ['mu_a_w', num2str(lamb), '_', num2str(frame_i), '.mat']);
    save(fname_mu_a, 'mu_a');
end
