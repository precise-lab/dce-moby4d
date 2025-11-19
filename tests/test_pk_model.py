import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../../hippylib/") )
sys.path.append("../")
import moby


time = np.arange(-60., 900.01, 0.05)
Ktrans  = 0.3/60 # 1/seconds
kep     = 0.75/60  # 1/seconds  
aif = moby.AIF()
pkModel = moby.PKModel(aif, time, Ktrans, kep)

aif_peak_idx = np.argmax(pkModel.ca_blood)
print("Peak AIF: ", pkModel.time[aif_peak_idx], " ", pkModel.ca_blood[aif_peak_idx])

plt.figure()
plt.plot(pkModel.time, pkModel.ca_blood, label="aif")
plt.plot(pkModel.time, pkModel.ca_viabletumor_perf, label="tumor")
plt.plot(pkModel.time, pkModel.ca_liver_perf, label="liver")
plt.savefig("pkModel.png")