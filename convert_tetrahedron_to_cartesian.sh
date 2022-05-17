#!/bin/bash

export HDF5_USE_FILE_LOCKING=FALSE

# Paths
export PATH=/software/paraview/bin:${PATH}

for frame_i in {1..1224}
do
  if (($frame_i % 10))
  then 
    echo moby_$(($frame_i % 10)).mat 
    pvpython convert_tetrahedron_to_cartesian.py \
      --phanfile /shared/planck/Phantom/Dynamic_MOBY/anatomical_structure_body_lesion/moby_$(($frame_i % 10)).mat \
      --phanvar phan \
      --tetfile /shared/planck/Phantom/Dynamic_MOBY/p0/p0_w800_${frame_i}.xdmf \
      --tetvar p0
  else
    echo moby_10.mat
    pvpython convert_tetrahedron_to_cartesian.py \
      --phanfile /shared/planck/Phantom/Dynamic_MOBY/anatomical_structure_body_lesion/moby_10.mat \
      --phanvar phan \
      --tetfile /shared/planck/Phantom/Dynamic_MOBY/p0/p0_w800_${frame_i}.xdmf \
      --tetvar p0
  fi
done
