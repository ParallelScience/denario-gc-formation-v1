# filename: codebase/step_2.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import interpolate
from astropy.cosmology import FlatLambdaCDM
import warnings

warnings.filterwarnings('ignore')

H0 = 67.4
Om0 = 0.315
Ob0 = 0.049
sigma8_target = 0.811
ns_spec = 0.965
h_hub = H0 / 100.0
Tcmb = 2.7255

cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0)

G_SI = 6.674e-11
Msun_kg = 1.989e30
Mpc_m = 3.0857e22
mu_mol = 0.6

rho_crit0_SI = 3.0 * (H0 * 1e3 / Mpc_m)**2 / (8.0 * np.pi * G_SI)
rho_crit0_Msun_Mpc3 = rho_crit0_SI * Mpc_m**3 / Msun_kg
rho_m0_Msun_Mpc3 = Om0 * rho_crit0_Msun_Mpc3
rho_m0_hMpc3 = rho_m0_Msun_Mpc3 / h_hub**3

def virial_temperature_to_mass_BL01(T_vir_K, z, Omega_m=Om0, h=h_hub):
    Delta_c = 18.0 * np.pi**2
    factor = 1.98e4 * (mu_mol / 0.6) * (Omega_m * Delta_c / (18.0 * np.pi**2))**(1.0 / 3.0) * ((1.0 + z) / 10.0)
    M_h_inv_Msun = (T_vir_K / factor)**(3.0 / 2.0)
    M_Msun = M_h_inv_Msun * 1e8 / h
    return M_Msun

if __name__ == '__main__':
    for z in [5, 10, 15, 20]:
        m_min = virial_temperature_to_mass_BL01(1e4, z)
        print('z=' + str(z) + ': M_min = ' + str(m_min) + ' M_sun')