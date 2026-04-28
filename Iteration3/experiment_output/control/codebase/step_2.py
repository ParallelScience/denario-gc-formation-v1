# filename: codebase/step_2.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import integrate, interpolate
from astropy.cosmology import FlatLambdaCDM
import os

H0 = 67.4
Om0 = 0.315
Ob0 = 0.049
sigma8_val = 0.811
ns_val = 0.965
h = H0 / 100.0
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0)

ST_a = 0.707
ST_p = 0.3
ST_A = 0.3222
mu = 0.59
T_vir_min = 1.0e4

rho_m_comoving = Om0 * 2.775e11 * h**2

def E_z(z):
    return np.sqrt(Om0 * (1.0 + z)**3 + (1.0 - Om0))

def delta_c_vir(z):
    x = Om0 * (1.0 + z)**3 / E_z(z)**2 - 1.0
    return 18.0 * np.pi**2 + 82.0 * x - 39.0 * x**2

def M_min_atomic_cooling(z):
    Delta_c = delta_c_vir(z)
    factor_delta = (Om0 / 0.3 * Delta_c / (18.0 * np.pi**2))**(1.0 / 3.0)
    factor_z = (1.0 + z) / 10.0
    T_coeff = 1.98e4 * (mu / 0.6) * factor_delta * factor_z
    M_h = (1.0e8 / h) * (T_vir_min / T_coeff)**(1.5)
    return M_h

def transfer_function_BBKS(k, Gamma):
    q = k / Gamma
    T = (np.log(1.0 + 2.34 * q) / (2.34 * q) *
         (1.0 + 3.89 * q + (16.1 * q)**2 + (5.46 * q)**3 + (6.71 * q)**4)**(-0.25))
    return T

def power_spectrum_shape(k, ns_v, Gamma):
    T = transfer_function_BBKS(k, Gamma)
    return k**ns_v * T**2

def compute_sigma_sq_R_unnorm(R, ns_v, Gamma):
    def integrand(lnk):
        k = np.exp(lnk)
        x = k * R
        W = 3.0 * (np.sin(x) - x * np.cos(x)) / x**3
        return power_spectrum_shape(k, ns_v, Gamma) * W**2 * k**3
    result, _ = integrate.quad(integrand, np.log(1e-4), np.log(1e4), limit=200)
    return result / (2.0 * np.pi**2)

def growth_factor(z):
    def integrand(a):
        z_a = 1.0 / a - 1.0
        return 1.0 / (a * E_z(z_a))**3
    norm, _ = integrate.quad(integrand, 1e-6, 1.0, limit=200)
    val, _ = integrate.quad(integrand, 1e-6, 1.0 / (1.0 + z), limit=200)
    return E_z(z) * val / (E_z(0.0) * norm)

if __name__ == '__main__':
    Gamma = Om0 * h * np.exp(-Ob0 * (1.0 + np.sqrt(2.0 * h) / Om0))
    sigma_sq_8_unnorm = compute_sigma_sq_R_unnorm(8.0, ns_val, Gamma)
    A_norm = sigma8_val**2 / sigma_sq_8_unnorm
    N_mass = 200
    M_vals = np.logspace(5, 13, N_mass)
    R_vals = (3.0 * M_vals / (4.0 * np.pi * rho_m_comoving))**(1.0 / 3.0)
    sigma_vals_z0 = np.array([np.sqrt(A_norm * compute_sigma_sq_R_unnorm(R_i, ns_val, Gamma)) for R_i in R_vals])
    log_sigma_interp = interpolate.interp1d(np.log(M_vals), np.log(sigma_vals_z0), kind='cubic', fill_value='extrapolate')
    z_arr = np.linspace(2.0, 20.0, 200)
    D_arr = np.array([growth_factor(z) for z in z_arr])
    delta_c_collapse = 1.686
    rate = np.zeros(len(z_arr))
    for i, z in enumerate(z_arr):
        M_min = M_min_atomic_cooling(z)
        D_z = D_arr[i]
        mask = M_vals >= M_min
        if mask.sum() < 3:
            rate[i] = 0.0
            continue
        M_sub = M_vals[mask]
        sigma_sub_z0 = np.exp(log_sigma_interp(np.log(M_sub)))
        sigma_sub = sigma_sub_z0 * D_z
        nu = delta_c_collapse / sigma_sub
        nu2 = ST_a * nu**2
        f_nu = (ST_A * np.sqrt(2.0 * ST_a / np.pi) * nu * (1.0 + nu2**(-ST_p)) * np.exp(-nu2 / 2.0))
        lnM_sub = np.log(M_sub)
        lnsigma_sub = np.log(sigma_sub)
        dlnsigma_dlnM = np.gradient(lnsigma_sub, lnM_sub)
        dndlnM = f_nu * (rho_m_comoving / M_sub) * np.abs(dlnsigma_dlnM)
        n_total = integrate.simpson(dndlnM, x=lnM_sub)
        dVdz_sr = cosmo.differential_comoving_volume(z).value
        dVdz_full = dVdz_sr * 4.0 * np.pi
        rate[i] = n_total * dVdz_full
    total_integral = integrate.simpson(rate, x=z_arr)
    pdf = rate / total_integral
    cdf = integrate.cumulative_trapezoid(pdf, z_arr, initial=0.0)
    np.savez('data/trenti_model.npz', z=z_arr, pdf=pdf, cdf=cdf)
    print('Peak redshift of theoretical GC formation rate PDF: ' + str(round(z_arr[np.argmax(pdf)], 3)))
    print('Saved theoretical model to data/trenti_model.npz')