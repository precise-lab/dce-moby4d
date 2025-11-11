import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.append("../")
import moby


time = np.arange(-30., 450.01, 0.05)
Ktrans  = 0.3/60 # 1/seconds
kep     = 0.75/60  # 1/seconds  
aif = moby.AIF()
pkModel = moby.PKModel(aif, time, Ktrans, kep)

plt.figure()
plt.plot(pkModel.time, pkModel.ca_blood, label="aif")
plt.plot(pkModel.time, pkModel.ca_viabletumor_perf, label="tumor")
plt.plot(pkModel.time, pkModel.ca_liver_perf, label="liver")
plt.savefig("pkModel.png")