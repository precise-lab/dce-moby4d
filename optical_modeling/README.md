## Dependencies:

`FEniCS 2019.1`: parallel finite element solver [https://fenicsproject.org/](https://fenicsproject.org/)

```bash
conda create -n fenics-2019.1 -c conda-forge fenics==2019.1.0 matplotlib scipy jupyter h5py hdf5storage meshio==5.3.4 pygalmesh==0.10.7 cgal
```

> Note: pygalmesh is bugged. See PR:

`hIPPYlib`: Large-Scale Inverse Problem in python (for inversion):

```bash
git clone https://github.com/hippylib/hippylib.git
```
