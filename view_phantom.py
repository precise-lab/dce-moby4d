import matplotlib.pyplot as plt
import scipy.io as io
import numpy as np

f = io.loadmat('/workspace/shared_data/Moby_multi_wave/phantom_800/moby000001.mat')['mu_sp'][20:-20,20:-20]
print(f.shape)
plt.imshow(f)
plt.savefig('anat_map.png')

plt.clf()

ca = io.loadmat('/workspace/shared_data/Moby_multi_wave/phantom_800/contrast_agent_curve.mat')['ca']
plt.plot(np.squeeze(ca))
plt.savefig('contract_curve.png')

