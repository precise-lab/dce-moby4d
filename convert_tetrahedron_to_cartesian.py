'''
AUTHOR: Seonyeong Park

CREATED:  10.15.2021

DESCRIPTION:
  Convert tetrahedron grid to Cartesian grid

'''
import os
import time
import meshio
import argparse
import numpy as np
import scipy.io as sio
from paraview.simple import *  # import the simple module from the paraview

if __name__ == '__main__':
  # Parse arguments
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--phanfile',
                      type = str,
                      help = 'Phantom .mat file name')
  parser.add_argument('--phanvar',
                      type = str,
                      help = 'Phantom variable name')
  parser.add_argument('--tetfile',
                      type = str,
                      help = 'Tetrahedron map .xdmf file name')
  parser.add_argument('--tetvar',
                      type = str,
                      help = 'Tetrahedron map variable name')
  args = parser.parse_args()

  # Measure time
  start = time.time()

  # Load 3D numerical phantom
  print('Loading 3D phantom .mat file...')
  data = sio.loadmat(args.phanfile)
  volume = data[args.phanvar]
  volume = volume[:,:,::-1] # THIS IS ONLY FOR MOBY P0!! (Z AXIS IS FLIPPED.)

  # Get sampling dimensions for Cartesian grid
  print('Calculating sampling dimensions for Cartesian grid...')
  Nx, Ny, Nz = volume.shape

  amax_x = np.amax(volume, axis=(1, 2))
  amax_y = np.amax(volume, axis=(0, 2))
  amax_z = np.amax(volume, axis=(0, 1))

  Nx_s = np.where(amax_x != 0)[0][0]
  Nx_e = np.where(amax_x != 0)[0][-1] + 1
  Ny_s = np.where(amax_y != 0)[0][0]
  Ny_e = np.where(amax_y != 0)[0][-1] + 1
  Nz_s = np.where(amax_z != 0)[0][0]
  Nz_e = np.where(amax_z != 0)[0][-1] + 1

  Nx_samp = Nx_e - Nx_s + 2
  Ny_samp = Ny_e - Ny_s + 2
  Nz_samp = Nz_e - Nz_s + 2

  # Load tetrahedron map (e.g., oxygen saturation disrtibution and initial pressure distribution)
  # create a new 'XDMF Reader'
  print('Loading tetrahedron map .xdmf file...')
  xdmf = XDMFReader(FileNames=args.tetfile)

  result = {}
  # Properties modified on xdmf
  xdmf.PointArrayStatus = args.tetvar
  xdmf.GridStatus = 'mesh'
  # Create a new 'Resample To Image'
  resampleToImage = ResampleToImage(Input=xdmf)
  # Properties modified on resampleToImage
  resampleToImage.SamplingDimensions = [Nx_samp, Ny_samp, Nz_samp]
  # Create a new 'Pass Arrays'
  passArrays = PassArrays(Input=resampleToImage)
  # passArrays.PointDataArrays = ['vtkValidPointMask', prop]
  # Properties modified on passArrays
  passArrays.PointDataArrays = args.tetvar
  # Save data as .vtk file
  filename_vtk = args.tetfile[:-5] + '.vtk'
  print('Saving ' + filename_vtk + '...')
  SaveData(filename_vtk, proxy=passArrays)

  # Load .vtk file
  print('Loading ' + filename_vtk + '...')
  mesh = meshio.read(filename_vtk)

  # Order correction for .mat file
  mesh_pd = mesh.point_data[args.tetvar].reshape(Nz_samp, Ny_samp, Nx_samp)
  mesh_pd = mesh_pd[::-1, :, :]
  mesh_pd = mesh_pd.transpose(2, 1, 0)

  # Crop domain to match with that of original phantom .mat file
  cartmap = np.zeros((Nx, Ny, Nz), dtype=np.float32)
  cartmap[Nx_s:Nx_e, Ny_s:Ny_e, Nz_s:Nz_e] = mesh_pd[1:-1, 1:-1, 1:-1]
  cartmap = cartmap[:,:,:-1]
  cartmap = cartmap[:,:,::-1] # THIS IS ONLY FOR MOBY P0!! (Z AXIS IS FLIPPED.)


  filename_cart_mat = args.tetfile[:-5] + '.mat'
  print('Saving ' + filename_cart_mat + '...')
  sio.savemat(filename_cart_mat, {args.tetvar: cartmap}, do_compression=True)

  # Delete vtk file if keep_vtk is False
#  print('Deleting ' + filename_vtk + '...')
#  os.system('rm ' + filename_vtk)

  print('  Elapsed time: %s sec' %(time.time() - start))

