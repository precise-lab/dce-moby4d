% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     May 3, 2022
%

function mu_sp = assign_mu_sp(properties, label_map, lamb)

% Reduced scattering coefficient mu_sp [1/mm]
mu_sp_table = properties.mu_s_ref.*(lamb./properties.wavelength_ref).^(-0.7).*(1 - properties.g);
mu_sp = mu_sp_table(label_map+1);
end