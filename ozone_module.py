"""
Plug-in ozone layer.  Give it ANY wavelength grid (nm); it returns the transmission T(lambda).

    from ozone_module import ozone_transmission
    T_o3 = ozone_transmission(lam_nm, temp=213, du=300, mu0=1.0)
    F_out = F_blackbody * T_o3 * T_your_N2          # layers simply multiply

Physics: T = exp(-N*sigma(lambda)/mu0), N = du*2.687e16 cm^-2, sigma in cm^2 (Serdyuchenko 2014).
Outside 213.33-1100 nm there is no data, so T = 1 there (NOT a claim that ozone is transparent).
temp must be one of the datasets: 193, 213 or 263 K.
"""
import os
import numpy as np
from ozone_layer import load_data, sigma_o3, tau, DU

_DATA = None

def ozone_transmission(lam_nm, temp=213, du=300.0, mu0=1.0, data_dir=None):
    global _DATA
    if _DATA is None:
        d = data_dir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        _DATA = load_data(d)[1]
    if temp not in _DATA:
        raise ValueError(f"temp must be one of {sorted(_DATA)} K")
    lam = np.asarray(lam_nm, float)
    return np.exp(-tau(du * DU, sigma_o3(lam, _DATA[temp])) / mu0)

if __name__ == "__main__":
    lam = np.array([200., 255.4, 300., 330., 450., 602., 800.])
    for T in (193, 213, 263):
        print(T, "K", np.array2string(ozone_transmission(lam, T), precision=3))
