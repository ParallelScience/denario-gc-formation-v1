# filename: codebase/step_1.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
from scipy import optimize
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import os

def get_cosmology():
    return FlatLambdaCDM(H0=67.4, Om0=0.315, Ob0=0.049)

def age_to_zform(z_obs, age_gyr, cosmo, z_min=0.001, z_max=100.0):
    t_obs = cosmo.lookback_time(z_obs).to(u.Gyr).value
    t_target = t_obs + age_gyr
    t_universe = cosmo.age(0).to(u.Gyr).value
    if t_target >= t_universe:
        return np.nan
    def func(z):
        return cosmo.lookback_time(z).to(u.Gyr).value - t_target
    try:
        return optimize.brentq(func, z_min, z_max, xtol=1e-6, rtol=1e-6)
    except ValueError:
        return np.nan

def compute_zform_with_errors(z_obs, age, age_err_lo, age_err_hi, cosmo):
    z_central = age_to_zform(z_obs, age, cosmo)
    z_lo_age = age_to_zform(z_obs, max(age - age_err_lo, 1e-6), cosmo)
    z_hi_age = age_to_zform(z_obs, age + age_err_hi, cosmo)
    if np.isnan(z_central):
        return np.nan, np.nan, np.nan
    err_lo = z_central - z_lo_age if not np.isnan(z_lo_age) else np.nan
    err_hi = z_hi_age - z_central if not np.isnan(z_hi_age) else np.nan
    return z_central, err_lo, err_hi

def tvir_to_mhalo(T_vir_K, z, cosmo, mu=0.59, Delta_c=200.0):
    G_cgs = 6.674e-8
    k_B_cgs = 1.3806e-16
    m_p_cgs = 1.6726e-24
    Msun_cgs = 1.989e33
    rho_crit_astropy = cosmo.critical_density(z).to(u.g / u.cm**3).value
    factor = (2.0 * k_B_cgs * T_vir_K) / (mu * m_p_cgs * G_cgs)
    M_halo_cgs = factor**(1.5) * (3.0 / (4.0 * np.pi * Delta_c * rho_crit_astropy))**(0.5)
    return M_halo_cgs / Msun_cgs

if __name__ == '__main__':
    cosmo = get_cosmology()
    data_dir = 'data/'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    print('Cosmology initialized. Ready for analysis.')