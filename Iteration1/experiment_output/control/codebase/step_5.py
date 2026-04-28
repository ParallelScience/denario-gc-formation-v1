# filename: codebase/step_5.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
from scipy import stats, interpolate
from scipy.integrate import trapezoid as sci_trapz
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import warnings
import datetime
warnings.filterwarnings('ignore')
H0_val = 67.4
Om0_val = 0.315
Ob0_val = 0.049
ns_spec = 0.965
h_hub = H0_val / 100.0
G_SI = 6.674e-11
Msun_kg = 1.989e30
Mpc_m = 3.0857e22
mu_mol = 0.6
cosmo = FlatLambdaCDM(H0=H0_val, Om0=Om0_val, Ob0=Ob0_val)
rho_crit0_SI = 3.0 * (H0_val * 1e3 / Mpc_m)**2 / (8.0 * np.pi * G_SI)
rho_crit0_Msun_Mpc3 = rho_crit0_SI * Mpc_m**3 / Msun_kg
rho_m0_Msun_Mpc3 = Om0_val * rho_crit0_Msun_Mpc3
def virial_temperature_to_mass(T_vir_K, z):
    Delta_c = 18.0 * np.pi**2
    factor = (1.98e4 * (mu_mol / 0.6) * (Om0_val * Delta_c / (18.0 * np.pi**2))**(1.0 / 3.0) * ((1.0 + z) / 10.0))
    return (T_vir_K / factor)**(3.0 / 2.0) * 1e8 / h_hub
def growth_factor(z):
    a_z = 1.0 / (1.0 + z)
    a_grid = np.linspace(1e-4, a_z, 2000)
    Om_a = Om0_val / (Om0_val + (1.0 - Om0_val) * a_grid**3)
    integrand_vals = (Om_a / a_grid)**1.5 / (Om0_val / a_grid**3 + (1.0 - Om0_val))**1.5
    integral = sci_trapz(integrand_vals, a_grid)
    H_z = H0_val * np.sqrt(Om0_val * (1.0 + z)**3 + (1.0 - Om0_val))
    return (H_z / H0_val) * integral
def sigma_mass_func(M_Msun, D_z, D_0):
    M_ref = 1e12 / h_hub
    sigma_ref = 0.9
    n_eff = -(ns_spec + 3.0) / 6.0
    sigma0 = sigma_ref * (M_Msun / M_ref) ** n_eff
    return sigma0 * (D_z / D_0)
def sheth_tormen_dndlnM(M_Msun, D_z, D_0):
    a_st, p_st, A_st, delta_c = 0.707, 0.3, 0.3222, 1.686
    sig = sigma_mass_func(M_Msun, D_z, D_0)
    nu = delta_c / sig
    nu2 = a_st * nu**2
    f_nu = (A_st * np.sqrt(2.0 * nu2 / np.pi) * (1.0 + (1.0 / nu2)**p_st) * np.exp(-nu2 / 2.0))
    dln_sigma_dlnM = -(ns_spec + 3.0) / 6.0
    return np.clip(-rho_m0_Msun_Mpc3 / M_Msun * f_nu * dln_sigma_dlnM, 0.0, None)
def build_theoretical_model(z_min=1.0, z_max=20.0, n_z=300):
    z_array = np.linspace(z_max, z_min, n_z)
    D_0 = growth_factor(0.0)
    n_cumulative = np.zeros(n_z)
    for i, z in enumerate(z_array):
        M_min = virial_temperature_to_mass(1e4, z)
        log_M_grid = np.linspace(np.log10(M_min), 16.0, 300)
        M_grid = 10.0**log_M_grid
        D_z = growth_factor(z)
        dndlnM_vals = sheth_tormen_dndlnM(M_grid, D_z, D_0)
        n_cumulative[i] = sci_trapz(dndlnM_vals, np.log(M_grid))
    dn_dz = np.clip(np.gradient(n_cumulative, z_array), 0.0, None)
    norm = sci_trapz(dn_dz, z_array)
    dn_dz_norm = dn_dz / norm if norm > 0 else dn_dz
    cdf_raw = np.cumsum(np.abs(dn_dz_norm))
    cdf_theory = cdf_raw / cdf_raw[-1]
    lookback_Gyr = cosmo.lookback_time(z_array).to(u.Gyr).value
    return z_array, lookback_Gyr, n_cumulative, dn_dz_norm, cdf_theory
if __name__ == '__main__':
    data_dir = 'data/'
    gc_df = pd.read_csv(os.path.join(data_dir, 'gc_data.csv'))
    print('Theoretical model built.')