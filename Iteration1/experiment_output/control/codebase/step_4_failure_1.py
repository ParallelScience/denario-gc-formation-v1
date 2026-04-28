# filename: codebase/step_4.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
from scipy import stats, integrate, interpolate
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import warnings
import os

warnings.filterwarnings('ignore')

H0 = 67.4
Om0 = 0.315
Ob0 = 0.049
sigma8_target = 0.811
ns_spec = 0.965
h_hub = H0 / 100.0
G_SI = 6.674e-11
Msun_kg = 1.989e30
Mpc_m = 3.0857e22
mu_mol = 0.6
kB_erg = 1.3806e-16
mp_g = 1.6726e-24
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0)
rho_crit0_SI = 3.0 * (H0 * 1e3 / Mpc_m)**2 / (8.0 * np.pi * G_SI)
rho_crit0_Msun_Mpc3 = rho_crit0_SI * Mpc_m**3 / Msun_kg
rho_m0_Msun_Mpc3 = Om0 * rho_crit0_Msun_Mpc3

def virial_temperature_to_mass(T_vir_K, z):
    Delta_c = 18.0 * np.pi**2
    factor = 1.98e4 * (mu_mol / 0.6) * (Om0 * Delta_c / (18.0 * np.pi**2))**(1.0 / 3.0) * ((1.0 + z) / 10.0)
    return (T_vir_K / factor)**(3.0 / 2.0) * 1e8 / h_hub

def sigma_mass(M_Msun, z=0):
    M_ref = 1e12 / h_hub
    sigma_ref = 0.9
    n_eff = -(ns_spec + 3.0) / 6.0
    sigma0 = sigma_ref * (M_Msun / M_ref) ** n_eff
    D_z = growth_factor(z)
    D_0 = growth_factor(0.0)
    return sigma0 * (D_z / D_0)

def growth_factor(z):
    def integrand(a):
        Om_a = Om0 / (Om0 + (1.0 - Om0) * a**3)
        return (Om_a / a)**1.5 / (Om0 / a**3 + (1.0 - Om0))**1.5
    a_z = 1.0 / (1.0 + z)
    a_grid = np.linspace(1e-4, a_z, 2000)
    integrand_vals = np.array([integrand(a) for a in a_grid])
    integral = np.trapz(integrand_vals, a_grid)
    H_z = H0 * np.sqrt(Om0 * (1.0 + z)**3 + (1.0 - Om0))
    return (H_z / H0) * integral

def sheth_tormen_dndlnM(M_Msun, z):
    a_st, p_st, A_st, delta_c = 0.707, 0.3, 0.3222, 1.686
    sig = sigma_mass(M_Msun, z)
    nu = delta_c / sig
    nu2 = a_st * nu**2
    f_nu = A_st * np.sqrt(2.0 * nu2 / np.pi) * (1.0 + (1.0 / nu2)**p_st) * np.exp(-nu2 / 2.0)
    dln_sigma_dlnM = -(ns_spec + 3.0) / 6.0
    return -rho_m0_Msun_Mpc3 / M_Msun * f_nu * dln_sigma_dlnM

def compute_cumulative_halo_density(z_array):
    n_cumulative = np.zeros(len(z_array))
    for i, z in enumerate(z_array):
        M_min = virial_temperature_to_mass(1e4, z)
        log_M_grid = np.linspace(np.log10(M_min), 16, 300)
        M_grid = 10.0**log_M_grid
        dndlnM_vals = np.clip(sheth_tormen_dndlnM(M_grid, z), 0.0, None)
        n_cumulative[i] = np.trapz(dndlnM_vals, np.log(M_grid))
    return n_cumulative

def compute_differential_rate(z_array, n_cumulative):
    dn_dz = np.clip(np.gradient(n_cumulative, z_array), 0.0, None)
    norm = np.trapz(dn_dz, z_array)
    return dn_dz / norm if norm > 0 else dn_dz

if __name__ == '__main__':
    data_dir = 'data/'
    df_gems = pd.read_csv(os.path.join(data_dir, 'gems_data.csv'))
    df_sparkler = pd.read_csv(os.path.join(data_dir, 'sparkler_data.csv'))
    z_array, _, _, _, cdf_theory = build_theoretical_model = (lambda z_min=1.0, z_max=20.0, n_z=200: (lambda z_arr: (z_arr, None, None, None, np.cumsum(np.abs(compute_differential_rate(z_arr, compute_cumulative_halo_density(z_arr)))) / np.cumsum(np.abs(compute_differential_rate(z_arr, compute_cumulative_halo_density(z_arr))))[-1]))(np.linspace(z_max, z_min, n_z)))()
    
    gems_z = df_gems['z_form'].values
    sparkler_z = df_sparkler['z_form'].values
    combined_z = np.concatenate([gems_z, sparkler_z])
    
    ad_stat, _, _ = stats.anderson_ksamp([gems_z, sparkler_z])
    
    z_theory_inc = z_array[::-1]
    cdf_theory_inc = cdf_theory[::-1]
    interp_cdf = interpolate.interp1d(z_theory_inc, cdf_theory_inc, bounds_error=False, fill_value=(0.0, 1.0))
    ks_stat, ks_p = stats.kstest(combined_z, interp_cdf)
    
    print('Anderson-Darling Test (GEMS vs Sparkler): Stat=' + str(ad_stat))
    print('One-sample KS Test (Combined vs Theory): Stat=' + str(ks_stat) + ', p=' + str(ks_p))
    
    for name, data in [('GEMS', df_gems), ('Sparkler', df_sparkler), ('Combined', pd.concat([df_gems, df_sparkler]))]:
        rho, p = stats.spearmanr(data['log10_mass'], data['z_form'])
        print(name + ' Spearman correlation (mass vs z_form): rho=' + str(rho) + ', p=' + str(p))