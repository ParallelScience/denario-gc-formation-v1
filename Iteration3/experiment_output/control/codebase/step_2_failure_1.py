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
sigma8 = 0.811
ns = 0.965
h = H0 / 100.0
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0)

ST_a = 0.707
ST_p = 0.3
ST_A = 0.3222
mu = 0.59
T_vir_min = 1.0e4

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
    T = np.log(1.0 + 2.34 * q) / (2.34 * q) * (1.0 + 3.89 * q + (16.1 * q)**2 + (5.46 * q)**3 + (6.71 * q)**4)**(-0.25)
    return T

def power_spectrum(k, ns_val, Gamma):
    T = transfer_function_BBKS(k, Gamma)
    return k**ns_val * T**2

def sigma_sq_R(R, A_norm, ns_val, Gamma):
    def integrand(lnk):
        k = np.exp(lnk)
        x = k * R
        W = 3.0 * (np.sin(x) - x * np.cos(x)) / x**3
        return power_spectrum(k, ns_val, Gamma) * W**2 * k**3
    result, _ = integrate.quad(integrand, np.log(1e-4), np.log(1e4), limit=200)
    return A_norm * result / (2.0 * np.pi**2)

def growth_factor(z):
    def integrand(a):
        z_a = 1.0 / a - 1.0
        return 1.0 / (a * E_z(z_a))**3
    norm, _ = integrate.quad(integrand, 0.0, 1.0, limit=100)
    val, _ = integrate.quad(integrand, 0.0, 1.0 / (1.0 + z), limit=100)
    return E_z(z) * val / (E_z(0.0) * norm)

if __name__ == '__main__':
    Gamma = Om0 * h * np.exp(-Ob0 * (1.0 + np.sqrt(2.0 * h) / Om0))
    sigma_sq_8_unnorm, _ = integrate.quad(lambda lnk: power_spectrum(np.exp(lnk), ns, Gamma) * (3.0 * (np.sin(np.exp(lnk)*8.0) - np.exp(lnk)*8.0 * np.cos(np.exp(lnk)*8.0)) / (np.exp(lnk)*8.0)**3)**2 * np.exp(lnk)**3, np.log(1e-4), np.log(1e4))
    A_norm = sigma8**2 / (sigma_sq_8_unnorm / (2.0 * np.pi**2))
    M_vals = np.logspace(5, 12, 100)
    sigma_vals = np.array([np.sqrt(sigma_sq_R((3.0 * M / (4.0 * np.pi * Om0 * 2.775e11 * h**2))**(1.0/3.0), A_norm, ns, Gamma)) for M in M_vals])
    sigma_interp = interpolate.interp1d(np.log(M_vals), np.log(sigma_vals), kind='cubic')
    z_arr = np.linspace(2, 20, 100)
    rate = []
    for z in z_arr:
        M_min = M_min_atomic_cooling(z)
        D_z = growth_factor(z)
        delta_c = 1.686
        mask = M_vals >= M_min
        if mask.sum() < 2:
            rate.append(0.0)
            continue
        M_sub = M_vals[mask]
        sigma_sub = np.exp(sigma_interp(np.log(M_sub))) * D_z
        nu = delta_c / sigma_sub
        f_nu = ST_A * np.sqrt(2.0 * ST_a / np.pi) * (1.0 + (ST_a * nu**2)**(-ST_p)) * np.sqrt(ST_a * nu**2) * np.exp(-ST_a * nu**2 / 2.0)
        dlnsigma_dlnM = np.gradient(np.log(sigma_sub), np.log(M_sub))
        dndlnM = f_nu * (Om0 * 2.775e11 * h**2 / M_sub) * np.abs(dlnsigma_dlnM)
        n_total = integrate.simps(dndlnM, np.log(M_sub))
        dVdz = cosmo.differential_comoving_volume(z).value * (h**3)
        rate.append(n_total * dVdz)
    rate = np.array(rate)
    pdf = rate / integrate.simps(rate, z_arr)
    cdf = integrate.cumtrapz(pdf, z_arr, initial=0)
    np.savez('data/trenti_model.npz', z=z_arr, pdf=pdf, cdf=cdf)
    print('Peak redshift:', z_arr[np.argmax(pdf)])