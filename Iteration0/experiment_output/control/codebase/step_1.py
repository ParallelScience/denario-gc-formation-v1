# filename: codebase/step_1.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.integrate import quad
from scipy.interpolate import interp1d
import os
import warnings

warnings.filterwarnings("ignore")

H0_val = 67.4
Om0_val = 0.315
Ob0_val = 0.049
OL0_val = 0.685
sigma8_val = 0.811
ns_val = 0.965

cosmo = FlatLambdaCDM(H0=H0_val, Om0=Om0_val, Ob0=Ob0_val)

G_SI = 6.674e-11
Msun_kg = 1.989e30
Mpc_m = 3.0857e22
kB = 1.381e-23
mp = 1.673e-27
mu = 0.6

ST_a = 0.707
ST_p = 0.3
ST_A = 0.3222

def bbks_transfer(k, cosmo_obj, ns):
    h = cosmo_obj.H0.value / 100.0
    Gamma = cosmo_obj.Om0 * h * np.exp(-cosmo_obj.Ob0 * (1.0 + np.sqrt(2.0 * h) / cosmo_obj.Om0))
    q = k / Gamma
    T = (np.log(1.0 + 2.34 * q) / (2.34 * q)) * (1.0 + 3.89 * q + (16.1 * q) ** 2 + (5.46 * q) ** 3 + (6.71 * q) ** 4) ** (-0.25)
    return k ** ns * T ** 2

def compute_sigma8_norm(cosmo_obj, ns):
    R8 = 8.0
    def integrand_sigma8(lnk):
        k = np.exp(lnk)
        x = k * R8
        W = 3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3
        P = bbks_transfer(k, cosmo_obj, ns)
        return k ** 3 * P * W ** 2 / (2.0 * np.pi ** 2) * k
    result, _ = quad(integrand_sigma8, np.log(1e-4), np.log(1e4), limit=200)
    return sigma8_val ** 2 / result

def sigma_M(M_msun, cosmo_obj, ns, A_norm):
    h = cosmo_obj.H0.value / 100.0
    rho_m_cgs = cosmo_obj.Om0 * cosmo_obj.critical_density(0).to(u.Msun / u.Mpc**3).value
    R_Mpc = (3.0 * M_msun / (4.0 * np.pi * rho_m_cgs)) ** (1.0 / 3.0)
    R_h = R_Mpc * h
    def integrand(lnk):
        k = np.exp(lnk)
        x = k * R_h
        W = 1.0 if x < 1e-3 else 3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3
        P = A_norm * bbks_transfer(k, cosmo_obj, ns)
        return k ** 3 * P * W ** 2 / (2.0 * np.pi ** 2) * k
    result, _ = quad(integrand, np.log(1e-4), np.log(1e4), limit=200)
    return np.sqrt(result)

def growth_factor(z, cosmo_obj):
    def integrand(a):
        z_a = 1.0 / a - 1.0
        H_a = cosmo_obj.H(z_a).value
        return 1.0 / (a * H_a) ** 3
    a_obs = 1.0 / (1.0 + z)
    norm, _ = quad(integrand, 0.0, 1.0, limit=200)
    val, _ = quad(integrand, 0.0, a_obs, limit=200)
    return (cosmo_obj.H(z).value * val) / (cosmo_obj.H(0).value * norm)

def delta_c_BN98(z, cosmo_obj):
    x = cosmo_obj.Om(z) - 1.0
    return 18.0 * np.pi ** 2 + 82.0 * x - 39.0 * x ** 2

def M_min_atomic_cooling(z, cosmo_obj):
    T_vir_min = 1.0e4
    Delta_c = delta_c_BN98(z, cosmo_obj)
    rho_crit_SI = cosmo_obj.critical_density(z).to(u.kg / u.m**3).value
    M_vir_kg = ((2.0 * kB * T_vir_min / (mu * mp)) ** 1.5) * ((3.0 / (4.0 * np.pi * Delta_c * rho_crit_SI)) ** 0.5) * (G_SI ** (-1.5))
    return M_vir_kg / Msun_kg

if __name__ == '__main__':
    A_norm = compute_sigma8_norm(cosmo, ns_val)
    z_grid = np.linspace(0, 30, 500)
    M_min_grid = np.array([M_min_atomic_cooling(z, cosmo) for z in z_grid])
    
    def get_n_cum(z):
        M_min = M_min_atomic_cooling(z, cosmo)
        rho_m = cosmo.Om0 * cosmo.critical_density(0).to(u.Msun / u.Mpc**3).value
        D_z = growth_factor(z, cosmo)
        delta_c = 1.686
        def integrand(lnM):
            M = np.exp(lnM)
            sigma = sigma_M(M, cosmo, ns_val, A_norm)
            nu = delta_c / (sigma * D_z)
            nu2 = ST_a * nu ** 2
            f_nu = ST_A * np.sqrt(2.0 * ST_a / np.pi) * nu * (1.0 + (1.0 / nu2) ** ST_p) * np.exp(-nu2 / 2.0)
            return f_nu * (rho_m / M)
        res, _ = quad(integrand, np.log(M_min), np.log(1e16), limit=100)
        return res

    n_cum_grid = np.array([get_n_cum(z) for z in z_grid])
    lookback_time = cosmo.lookback_time(z_grid).value
    dndt = np.gradient(n_cum_grid, lookback_time)
    norm_rate = dndt / np.max(dndt)
    
    if not os.path.exists("data"):
        os.makedirs("data")
    np.savez("data/theoretical_model.npz", z=z_grid, n_cum=n_cum_grid, M_min=M_min_grid, rate=norm_rate)
    print("Theoretical model saved to data/theoretical_model.npz")