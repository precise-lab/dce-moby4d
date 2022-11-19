% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     May 3, 2022
%

function mu_a = assign_mu_a(properties, label_map, lamb, ca, c_perf)


f_b = properties.f_b_oxy + properties.f_b_deoxy;

f_ca = f_b*ca;
chi_lesion = strcmp(properties.label2str(:,1), 'lesn');
f_ca(chi_lesion) = f_ca(chi_lesion) + (1. - f_b(chi_lesion))*c_perf;

mu_a_table = properties.f_b_oxy*properties.mu_a_b_oxy(properties.wavelength == lamb) + ...
    properties.f_b_deoxy*properties.mu_a_b_deoxy(properties.wavelength == lamb) + ...
    properties.f_w*properties.mu_a_w(properties.wavelength == lamb) + ...
    f_ca*properties.mu_a_ca(properties.wavelength == lamb);

mu_a = mu_a_table(label_map+1);

end