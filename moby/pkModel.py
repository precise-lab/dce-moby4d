import numpy as np

class AIF:
    def __init__(self, scaling: float = 5):
        self.A = scaling * 120 * 1e-6  # M/ min^Ap
        self.Ap = 3
        self.B = 4.34/60    # 1 / seconds
        self.C = scaling * 0.8*1e-6 # M
        self.D = 1/60   # 1/ seconds
        self.E = 0.07/60 # 1/ seconds
        self.T0 = 0. #seconds

    def eval(self, time):
        out = np.zeros_like(time)
        t = time[time > self.T0]
        out[time > self.T0] = self.A*((t/60.)**self.Ap)*np.exp(-self.B*t) + self.C*(1 - np.exp(-self.D*t))*np.exp(-self.E*t)
        return out

class PKModel:
    def __init__(self, aif: AIF, time: np.array, Ktrans, kep):
        self.aif  = aif
        self.Ktrans = Ktrans
        self.kep = kep
        dt = time[1] - time[0]
        self.time = time
        self.ca_blood = self.aif.eval(time)
        aif_peak_index = np.argmax(self.ca_blood)
        self.aif_peak_time = self.time[aif_peak_index]
        t = time[time > aif.T0]
        self.ca_viabletumor_perf = np.zeros_like(self.ca_blood)
        self.ca_viabletumor_perf[time > aif.T0] = (self.Ktrans*np.convolve(np.exp(-self.kep*t), self.ca_blood[time > aif.T0])*dt)[:t.shape[0]]
        self.ca_core_perf        = np.zeros_like(self.ca_blood)
        self.ca_core_perf[time > aif.T0]        = (self.Ktrans*np.convolve(np.exp(-self.kep*t), self.ca_blood[time > aif.T0])*dt)[:t.shape[0]]

        self.ca_liver_perf                          = np.zeros_like(self.ca_blood)
        self.ca_liver_perf[time>self.aif_peak_time] = np.maximum(2.*(aif.C - self.ca_blood[time>self.aif_peak_time]), 0)
        

    def __call__(self, t: float):
        return [np.interp(t, self.time, self.ca_blood),
                np.interp(t, self.time, self.ca_viabletumor_perf),
                np.interp(t, self.time, self.ca_core_perf),
                np.interp(t, self.time, self.ca_liver_perf)]


