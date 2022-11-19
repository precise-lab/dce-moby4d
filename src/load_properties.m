function properties = load_properties(prefix,f_ICG_in_blood)

% Load functional properties
load(fullfile(prefix, 'func_prop.mat'));

properties.f_b_oxy     = s.*f_b;
properties.f_b_deoxy   = (1.-s).*f_b;
properties.f_w = f_w;

% Load tissue scattering properties
load(fullfile(prefix, 'opt_prop.mat'))
properties.mu_s_ref = mu_s_ref;
properties.g         = g;
properties.wavelength_ref = wavelength_ref;


% Load constants
load(fullfile(prefix, 'constants.mat'), ...
    'wavelength', 'e_hbo2', 'e_hb', 'mu_a_w');

properties.wavelength = wavelength;
properties.mu_a_b_oxy = log(10)*c_thb_b*e_hbo2;
properties.mu_a_b_deoxy = log(10)*c_thb_b*e_hb;
properties.mu_a_w = mu_a_w;


% Load molar extinction coefficient of ICG
load(fullfile(prefix, 'constant_icg.mat'), 'e_icg');
properties.mu_a_ca = f_ICG_in_blood*log(10)*e_icg(:,1)*1e-3;

load(fullfile(prefix, 'label.mat'));
properties.label2str = label;


end