# Dynamic PACT MOBY phantom
This project is to generate frames of MOBY mouse and spherical lesion phantoms with motion (heart beating and/or respiratory motion) and assign functional and optical properties to each tissue type.


## 1. MOBY body and lesion phantoms
### Configuration
The phantoms are generated with the following configuration.
- Voxel size: 0.15 mm
- Dimension size: 248 x 248 x 725
- Heart beating cycle: 0.175 s
- Respiratory cycle: 1.025 s
- Time per frame: 0.05 s
- Number of frames: 287 

All parameters used to generate MOBY body and lesion phantoms are in `anatomical_structure_body/moby.par` and `anatomical_structure_lesion/lesion.par`, respectively.

### Generate body phantom
The following command line generates `anatomical_structure_body/moby_act_{FRAME_1_to_287}.bin` files.
```bash
./moby_64bit anatomical_structure_body/moby.par anatomical_structure_body/moby
```

### Generate spherical lesion phantom
The following command line generates `anatomical_structure_lesion/lesion_act_{FRAME_1_to_287}.bin` files
```bash
./moby_64bit anatomical_structure_lesion/lesion.par anatomical_structure_lesion/lesion
```

### Data location
Generated phantoms can be found at https://utexas.app.box.com/folder/167127747374?s=z3d9r7fkdiypa7j6ignphnz21rc1jyjk

## 2. Dynamic optical phantoms for DCE-PACT 
### Constant and property tables for optical property assignment
Before running codes to build the optical phanton, run `init_properties.m`.

This will create the following files:

- `properties/label.mat`: 147 x 3 cell array `label` (Column 1: MOBY tissue name; Column 2: Mouse atlas tissue name; Column 3: integer tissue label)
- `properties/func_prop.mat`: total hemoglobin concentration in blood `c_thb_b` (unit: mol/L, M) and 147 x 1 arrays of oxygen saturation `s`, blood volume fraction `f_b`, and water volume fraction `f_w`
- `properties/opt_prop.mat`: 147 x 1 arrays of reference wavelength `wavelength_ref` (unit: nm), scattering coefficient at the reference wavelength `mu_s_ref` (unit: 1/mm), and scattering anisotropy `g`
- `properties/constants.mat`: 571 x 1 array of wavelength `wavelength` (unit: nm), 571 x 1 arrays of molar extinction coefficient of deoxy- `e_hb` and oxyhemoglobin `e_hbo2` (unit: 1/(mmM)) and optical absorption coefficient of water `mu_a_w` (unit: 1/mm). Other constants are not used here.
- `properties/constant_icg.mat`: 571 x 1 array of wavelength `wavelength` (unit: nm), 571 x 4 array of molar extinction coefficient of ICG `e_icg` (unit: 1/(mmM); Column 1: concentration `c` of 6.5 μm; Column 2: `c` of 65 μm; Column 3: `c` of 650 μm; Column 4: `c` of 1290 μm)


### Build the optical phantom for the DCE experiment
Run `driver.m`, this will generate the `config file` for use in step 3, and 4D images of tissue label maps, optical absorption and reduced scattering distribution.

## 3. Simulation of initial pressure distribution

### Generation of unstructured meshes for each anatomical frame
```
python optical_modeling/generate_meshes.py --config <PATH_TO_CONFIG_FILE>
``` 

### Simulation of the initial pressure distribution
```
mpirun -n 16 pythonoptical_modeling/generate_dynamic_p0.py --config <PATH_TO_CONFIG_FILE>
``` 

## Dependencies:

`FEniCS 2019.1`: parallel finite element solver [https://fenicsproject.org/](https://fenicsproject.org/)

```bash
conda create -n fenics-2019.1 -c conda-forge fenics==2019.1.0 matplotlib scipy jupyter h5py hdf5storage meshio==5.3.4 pygalmesh==0.10.7 cgal
```

> Note: pygalmesh is bugged. See PR: https://github.com/meshpro/pygalmesh/pull/204

`hIPPYlib`: Large-Scale Inverse Problem in python (for inversion):

```bash
git clone https://github.com/hippylib/hippylib.git
```