# Dynamic PACT MOBY phantom
This project is to generate frames of MOBY mouse and spherical lesion phantoms with motion (heart beating and/or respiratory motion) and assign functional and optical properties to each tissue type.


## MOBY body and lesion phatnoms
### Configuration
The phantoms are generated with the following configuration.
- Voxel size: 0.1 mm
- Dimension size: 372 x 372 x 1088
- Heart beating cycle: 1 s
- Respiratory cycle: 5 s
- Time per frame: 0.5 s
- Number of frames: 10 (1 respiratory cycle)

All parameters used to generate MOBY body and lesion phantoms are in `anatomical_structure_body/moby.par` and `anatomical_structure_lesion/lesion.par`, respectively.

### Generate body phantom
The following command line generates `anatomical_structure_body/moby_act_{FRAME_1_to_10}.bin` files (every 10 frames repeat due to a respiratory cycle of 5 s).
```bash
./moby_64bit anatomical_structure_body/moby.par anatomical_structure_body/moby
```

### Generate spherical lesion phantom
The following command line generates `anatomical_structure_lesion/lesion_act_{FRAME_1_to_10}.bin` files (every 10 frames repeat due to a respiratory cycle of 5 s).
```bash
./moby_64bit anatomical_structure_lesion/lesion.par anatomical_structure_lesion/lesion
```

### Insert the lesion to the body
The following MATLAB command generates `anatomical_structure_body_lesion/moby_{FRAME_1_to_10}.mat` files (every 10 frames repeat due to a respiratory cycle of 5 s).
```
run('insert_lesion_to_body.m');
```


## Constant and property tables for optical property assignment
Before running codes to assign optical absorption coefficient `mu_a` and reduced scattering coefficient `mu_sp`,
- Run `properties/save_label.m` to create the `properties/label.mat` file;
- Run `properties/save_func_prop.m` to create the `properties/func_prop.mat` file;
- Run `properties/save_opt_prop.m` to create the `properties/opt_prop.mat` file;
- Run `properties/constant_icg.m` to create the `properties/constant_icg` file;
- Run `properties/contrast_agent_curve.m` to create the `properties/contrast_agent_curve.mat` file.

### Constant and property table files
- `properties/label.mat`: 147 x 3 cell array `label` (Column 1: MOBY tissue name; Column 2: Mouse atlas tissue name; Column 3: integer tissue label)
- `properties/func_prop.mat`: total hemoglobin concentration in blood `c_thb_b` (unit: mol/L, M) and 147 x 1 arrays of oxygen saturation `s`, blood volume fraction `f_b`, and water volume fraction `f_w`
- `properties/opt_prop.mat`: 147 x 1 arrays of reference wavelength `wavelength_ref` (unit: nm), scattering coefficient at the reference wavelength `mu_s_ref` (unit: 1/mm), and scattering anisotropy `g`
- `properties/constants.mat`: 571 x 1 array of wavelength `wavelength` (unit: nm), 571 x 1 arrays of molar extinction coefficient of deoxy- `e_hb` and oxyhemoglobin `e_hbo2` (unit: 1/(mmM)) and optical absorption coefficient of water `mu_a_w` (unit: 1/mm). Other constants are not used here.
- `properties/constant_icg.mat`: 571 x 1 array of wavelength `wavelength` (unit: nm), 571 x 4 array of molar extinction coefficient of ICG `e_icg` (unit: 1/(mmM); Column 1: concentration `c` of 6.5 μm; Column 2: `c` of 65 μm; Column 3: `c` of 650 μm; Column 4: `c` of 1290 μm)
- `properties/contrast_agent_curve.mat`: 1 x 1224 arrays of time `t` (unit: s), concentration of contrast agent in plasma `cp`, and concentration of contrast agent in tissue of interest `ctoi`. Other parameters used to compute `cp` and `ctoi` (`Ktrans` and `kep` in two-compartment model, `A`, `B`, `C`, `D`, and `E` in arterial input function, etc.) are included.


## Optical property assignment
### Optical absorption coefficient
The following MATLAB command generates `mu_a/mu_a_{FRAME_1_to 1224}.mat` files (`mu_a` map varies over time regardless of a respiratory cycle of 5 s).
```
run('assign_mu_a.m');
```

### Reduced scattering coefficient
The following MATLAB command generates `mu_sp/mu_sp_{FRAME_1_to 10}.mat` files (every 10 frames repeat due to a respiratory cycle of 5 s).
```
run('assign_mu_sp.m');
```