import dolfin as dl
import ufl
import numpy as np

import argparse

from timeit import default_timer as timer

import sys
import os

sys.path.append("../")
import moby


sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib/") )
import hippylib as hp

class Material:
    def __init__(self, label, c_thb = 0., sat = 1., c_w = 0., c_f = 0., c_m = 0., mu_ps0 = 1.92, g = 0.95):
        self.label  = label
        self._c_thb = c_thb
        self._sat   = sat
        self._c_w   = c_w
        self._c_f   = c_f
        self._c_m   = c_m
        self._mu_ps0 = mu_ps0
        self._g     = g
        
    def to_dict(self):
        out = {}
        out["label"]  = self.label
        out["c_thb"]  = self._c_thb
        out["sat"]    = self._sat
        out["c_w"]    = self._c_w
        out["c_f"]    = self._c_f
        out["c_m"]    = self.c_m
        out["mu_ps0"] = self._mu_ps0
        out["g"]      = self._g
        
        return out
        
        
    
    def as_vector(self, vars=None):
        if vars is None:
            return dl.as_vector([self.c_thb*(dl.Constant(1.)-self.sat),
                             self.c_thb*self.sat,
                             self.c_w,
                             self.c_f,
                             self.c_m,
                             self.mu_ps0])
        else:
            tmp = []
            for var in vars:
                if var == Chromophores.hb:
                    tmp.append(self.c_thb*(dl.Constant(1.)-self.sat))
                elif var == Chromophores.hbO2:
                    tmp.append(self.c_thb*self.sat)
                elif var == Chromophores.water:
                    tmp.append(self.c_w)
                elif var == Chromophores.fat:
                    tmp.append(self.c_f)
                elif var == Chromophores.melanin:
                    tmp.append(self.c_m)
                elif var == "mu_ps0":
                    tmp.append(self.mu_ps0)
                elif var == "mu_s0":
                    tmp.append(self.mu_s0)
                elif var == "g":
                    tmp.append(self.g)
                elif var == "s":
                    tmp.append(self.sat)
                elif var == "thb":
                    tmp.append(self.c_thb)
                else:
                    raise ValueError(var)
                
            return dl.as_vector(tmp)

        
        
    @property
    def c_thb(self):
        return dl.Constant(self._c_thb)
    
    @property
    def sat(self):
        try:
            return dl.Constant(self._sat)
        except:
            return self._sat
        
    
    @property
    def c_w(self):
        return dl.Constant(self._c_w)
    
    @property
    def c_f(self):
        return dl.Constant(self._c_f)
    
    @property
    def c_m(self):
        return dl.Constant(self._c_m)
    
    @property
    def mu_ps0(self):
        return dl.Constant(self._mu_ps0)
    
    @property
    def mu_s0(self):
        return dl.Constant(self._mu_ps0/(1.-self._g))
    
    @property
    def g(self):
        return dl.Constant(self._g)

def define_materials(Labels):
    materials = {}

    
    materials[Labels['artery']]    = Material(label = Labels['artery'],
                                           c_thb = .99,
                                           sat   = .98,
                                           c_w   = 0.0,
                                           mu_ps0 = 1.3175715,
                                           g     = 0.9826)
    
    materials[Labels['tumor']]    = Material(label = Labels['tumor'],
                                           c_thb = .05,
                                           sat   = .371,
                                           c_w   = 0.0,
                                           mu_ps0 = 1.3175715,
                                           g     = 0.9826)

    materials[Labels['tumor_core']]    = Material(label = Labels['tumor_core'],
                                           c_thb = .0,
                                           sat   = .0,
                                           c_w   = 0.0,
                                           mu_ps0 = 1.3175715,
                                           g     = 0.9826)

    materials[Labels['vein']]      = Material(label = Labels['vein'],
                                           c_thb = .99,
                                           sat   = .7,
                                           c_w   = 0.0,
                                           mu_ps0 = 1.278304,
                                           g     = 0.9828)
    
    return materials


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Compute SO2', fromfile_prefix_chars='@')
    parser.add_argument('-f', '--fname', default = "/workspace/shared_data/Moby_multi_wave/mesh/")
    parser.add_argument('-o', '--output', default = "/workspace/shared_data/Moby_multi_wave/so2/")
    parser.add_argument('--start_frame',
                        default=0,
                        type=int,
                        help = "First frame index")
    parser.add_argument('--end_frame',
                        default=200,
                        type=int,
                        help = "End frame index")

    cl_args = parser.parse_args()

    comm = dl.MPI.comm_world
    tissueComposition = moby.TissueComposition.create()
    Labels = tissueComposition.tissue2label
    materials = define_materials(Labels)


    for i in range(cl_args.start_frame, cl_args.end_frame):
        print(i)
        mesh = dl.Mesh(dl.MPI.comm_world)
        if i > 0:
            with dl.XDMFFile(cl_args.fname + f"moby_mesh{i}.xdmf") as fid:
                fid.read(mesh)
                geo_dim = mesh.geometry().dim()
                c_labels = dl.MeshFunction('size_t', mesh, geo_dim)
                fid.read(c_labels, "c_labels")
        else:
            with dl.XDMFFile("/workspace/shared_data/Moby_multi_wave/mesh/moby_mesh_structured.xdmf") as fid:
                fid.read(mesh)
                geo_dim = mesh.geometry().dim()
                c_labels = dl.MeshFunction('size_t', mesh, geo_dim)
                fid.read(c_labels, "c_labels")
        dx = dl.Measure("dx", subdomain_data=c_labels)

        Vh = dl.FunctionSpace(mesh, 'CG', 1)  
        #Vh = dl.FunctionSpace(mesh, "DG", 0)
        sO2diff = dl.Constant(1e-6)

        uh, vh = dl.TrialFunction(Vh), dl.TestFunction(Vh)
        varf = sO2diff* dl.inner(dl.grad(uh), dl.grad(vh))*dl.dx \
          + uh*vh*dx(Labels['artery']) + uh*vh*dx(Labels['vein']) + uh*vh*dx(Labels['tumor']) + + uh*vh*dx(Labels['tumor_core'])
        rhs  = materials[Labels['artery']].sat*vh*dx(Labels['artery']) + \
           materials[Labels['vein']].sat*vh*dx(Labels['vein']) + \
           materials[Labels['tumor']].sat*vh*dx(Labels['tumor']) + \
           materials[Labels['tumor_core']].sat*vh*dx(Labels['tumor_core'])
     
        u = dl.Function(Vh, name='Saturation')
        
        A, b = dl.assemble_system(varf, rhs, [])
        
        dl.solve(A, u.vector(), b, 'cg', 'hypre_amg')

        with dl.XDMFFile(cl_args.output+f"so2_{i}.xdmf") as out:
            out.write_checkpoint(u, "so2", append=False)
            

