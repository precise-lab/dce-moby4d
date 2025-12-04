The current pipeline is contained in the director pipeline/

There is a shell script 
```bash
./run_pipeline.sh
```
that will run the full pipeline begginning with the following

- 'compute_displacement_field.py': computes the field of displacement vectors with inputes  from '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_vectors_converted_to_numpy_arrays/' and outputs the displacement field to  '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/motion_array/
- 'combine_phantom_elements.py': combines phantom elements for a single frame taking in base anatomy from '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/moby_anatomy/moby_act_1.bin', vein from '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/onlyvein_r15.DAT', artery from '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/onlyartery_r15.DAT', and tumor from '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/tumor_phantom.mat' then outputting single anatamy frame to '/workspace/shared_data/Moby_multi_wave/Refik_Mouse/phantom_anatomy_0.npy'
- 'generate_mesh_structured.py': generates structured mesh from single anatomy frame and saves it to '/workspace/shared_data/Moby_multi_wave/mesh/moby_mesh_structured.xdmf'
- 'move_mesh.py': moves mesh based on displacement field and saves dynamic mesh to '/workspace/shared_data/Moby_multi_wave/mesh/'
- 'compute_so2.py': computes oxygen saturation from mesh and saves it to '/workspace/shared_data/Moby_multi_wave/so2/'
- 'compute_mask.py': computes mask for arterys, tumor, and necrotic core then saves it to '/workspace/shared_data/Moby_multi_wave/masks/'

The pipeline then runs several scripts of the form 'run_...sh' which computes the optical forward model at a specificed wavelength then saves the resulting induced pressure via the file 

- 'combined_optical_forward.py': this comprises the majority of the compute time and phantom generation pipeline

Once the optical forwad is computed 

-'static_mask.py': computes a static mask for perfusion analysis based on the dynamic masks

To generate acoustic measurements utilize the directory 'acoustic_modeling/', with particular files of not being

- 'generate_measurements.py': which generates acoustic measurements at the subssampled resolution
- 'subsample_p0.py': computes a subsampled p0 for error calculation
