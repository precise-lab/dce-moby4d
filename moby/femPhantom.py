import dolfin as dl
import ufl
import csv
from .tissue_composition import TissueComposition
from .pkModel import PKModel
from .chromophore_decomposition import ChromophoreDecomposition, Chromophore

class FEMPhantom:
    def __init__(self, Vh_DG0, Vh_CG1, cell_labels,
                       tissue_composition: TissueComposition,
                       chromophore_decomposition: ChromophoreDecomposition,
                       pk_model: PKModel):
        self.Vh_DG0 = Vh_DG0
        self.Vh_CG1 = Vh_CG1
        self.cell_labels = cell_labels
        self.mesh = Vh_DG0.mesh()
        self.tissueComposition = tissue_composition
        self.chromophoreDecomposition = chromophore_decomposition
        self.pkModel = pk_model

        #Used chromophores
        self.chromophore_list = [Chromophore.HB, Chromophore.HBO2, Chromophore.WATER, Chromophore.CA]
        
        #define dx
        self.dx = dl.Measure("dx", subdomain_data=self.cell_labels, domain = self.mesh)

        #define vols
        labels = self.tissueComposition.tissue2label
        self.vols = [dl.assemble(dl.Constant(1.)*self.dx(labels[key])) for key in labels.keys()]

        #define mass matrix DG0
        trial_DG, test_DG = ufl.TrialFunction(self.Vh_DG0), ufl.TestFunction(self.Vh_DG0)
        varf = ufl.inner(trial_DG, test_DG)*self.dx
        self.M_DG0 = dl.assemble(varf)
        

    def compute_oxygen_saturation(self, tissue_oxygen_saturation: dict):

        labels = self.tissueComposition.tissue2label
        sO2diff = dl.Constant(1e-6)
        uh, vh = ufl.TrialFunction(self.Vh_P1), ufl.TestFunction(self.Vh_P1)
        varf = sO2diff* ufl.inner(ufl.grad(uh), ufl.grad(vh))*self.dx
        rhs = dl.Constant(0.)*vh*self.dx

        for tissue, sat in tissue_oxygen_saturation.items():
            varf += uh*vh*self.dx(labels[tissue])
            rhs  += dl.Constant(sat)*vh*self.dx(labels[tissue])

        out = dl.Function(self.Vh_CG1, name='Saturation')
        
        A, b = dl.assemble_system(varf, rhs, [])
        
        dl.solve(A, out.vector(), b, 'cg', 'hypre_amg')

        return out

    def compute_mu_sp(self, wavelength):
        labels = self.tissueComposition.tissue2label
        test = dl.TestFunction(self.Vh_DG0)
        rhs = dl.Constant(0.)*test*self.dx
        for key, vol in zip(labels.keys(), self.vols):
            if vol > 0:
                rhs += dl.Constant(self.tissue_composition.getReducedScattering(labels[key], wavelength)*test*self.dx(labels[key]))
        b = dl.assemble_system(rhs)
        out = dl.Function(self.Vh_DG, name = 'mu_sp')
        dl.solve(self.M_DG0, out.vector(), b, 'cg', 'jacobi')

        return out
    
    def compute_mu_a(self, wavelength, time, sO2):
        labels = self.tissueComposition.tissue2label
        test = dl.TestFunction(self.Vh_DG0)
        rhs = dl.Constant(0.)*test*self.dx


        aif, c_perf = self.pkModel.eval(time)
        c_thb_b = self.tissueComposition.c_thb_b
        c_CA_b  = self.pkModel.c_CA_b

        pure_mu_a = self.chromophoreDecomposition.get_pure_mu_a(wavelength, c_thb_b, c_CA_b, self.chromophore_list)

        for key, vol in zip(labels.keys(), self.vols):
            if vol > 0:
                f_b = self.tissueComposition.getVolumeFraction(labels[key], "blood")
                f_w = self.tissueComposition.getVolumeFraction(labels[key], "water")
                f_ca = f_b*aif
                if key == 'tumor':
                    f_ca +=  (1. - f_b)*c_perf
                if key == 'tumor_core':
                    f_ca +=  0.5*(1. - f_b)*c_perf
                        
                s_coeff = (pure_mu_a[Chromophore.HBO2]-pure_mu_a[Chromophore.HB])*f_b
                c_coeff = pure_mu_a[Chromophore.HB]*f_b \
                          + pure_mu_a[Chromophore.WATER]*f_w \
                          + pure_mu_a[Chromophore.CA]*f_ca
                rhs += (dl.Constant(s_coeff)*sO2 + dl.Constant(c_coeff))*test*self.dx(labels[key])
                
        b = dl.assemble_system(rhs)
        out = dl.Function(self.Vh_DG0, name = 'mu_a')
        dl.solve(self.M_DG0, out.vector(), b, 'cg', 'jacobi')

    def dumpTissueVolumes(self, fname):
        labels = self.tissueComposition.tissue2label
        fieldnames = ["tissue", "label", "volume"]
        with open(fname, "w") as fid:
            writer = csv.DictWriter(fid, fieldnames=fieldnames)
            writer.writeheader()
            rows = []
            for key, vol in zip(labels.keys(), self.vols):
                rows.append({"tissue": key, "label": labels[key], "volume": vol})
            writer.writerows(rows)
            
