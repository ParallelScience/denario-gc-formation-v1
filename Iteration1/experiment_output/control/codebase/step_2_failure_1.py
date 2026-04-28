# filename: codebase/step_2.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import integrate, interpolate
from astropy.cosmology import FlatLambdaCDM
import warnings

warnings.filterwarnings('ignore')

H0 = 67.4
Om0 = 0.315
Ob0 = 0.049
sigma8 = 0.811
ns = 0.965
h = H0 / 100.0
Tcmb = 2.7255
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0)
G_SI = 6.674e-11
Msun_kg = 1.989e30
Mpc_m = 3.0857e22
kB = 1.381e-23
mp = 1.673e-27
mu = 0.59
rho_crit0_SI = 3.0 * (H0 * 1e3 / Mpc_m)**2 / (8.0 * np.pi * G_SI)
rho_crit0_Msun_Mpc3 = rho_crit0_SI * Mpc_m**3 / Msun_kg
rho_m0_Msun_Mpc3 = Om0 * rho_crit0_Msun_Mpc3
delta_c = 1.686
ST_a = 0.707
ST_p = 0.3
ST_A = 0.3222

def eisenstein_hu_transfer(k_hMpc, Omega_m, Omega_b, h, Tcmb):
    k = k_hMpc * h
    ombh2 = Omega_b * h**2
    ommh2 = Omega_m * h**2
    Theta = Tcmb / 2.7
    z_eq = 2.5e4 * ommh2 * Theta**(-4)
    k_eq = 7.46e-2 * ommh2 * Theta**(-2)
    b1 = 0.313 * ommh2**(-0.419) * (1.0 + 0.607 * ommh2**0.674)
    b2 = 0.238 * ommh2**0.223
    z_d = 1291.0 * ommh2**0.251 / (1.0 + 0.659 * ommh2**0.828) * (1.0 + b1 * ombh2**b2)
    R_eq = 31.5e3 * ombh2 * Theta**(-4) * (1000.0 / z_eq)
    R_d = 31.5e3 * ombh2 * Theta**(-4) * (1000.0 / z_d)
    s = 2.0 / (3.0 * k_eq) * np.sqrt(6.0 / R_eq) * np.log((np.sqrt(1.0 + R_d) + np.sqrt(R_d + R_eq)) / (1.0 + np.sqrt(R_eq)))
    k_silk = 1.6 * (ombh2**0.52) * (ommh2**0.38) * (1.0 / (1.0 + (5.2 * ommh2)**(-1.36)))**(-1)
    a1 = (46.9 * ommh2)**0.670 * (1.0 + (32.1 * ommh2)**(-0.532))
    a2 = (12.0 * ommh2)**0.424 * (1.0 + (45.0 * ommh2)**(-0.582))
    alpha_c = a1**(-Omega_b / Omega_m) * a2**(-(Omega_b / Omega_m)**3)
    bb1 = 0.944 / (1.0 + (458.0 * ommh2)**(-0.708))
    bb2 = (0.395 * ommh2)**(-0.0266)
    beta_c = 1.0 / (1.0 + bb1 * ((Omega_b / Omega_m)**bb2 - 1.0))
    def T_tilde(k_arr, alpha, beta):
        q = k_arr / (13.41 * k_eq)
        C = 14.2 / alpha + 386.0 / (1.0 + 69.9 * q**1.08)
        T0 = np.log(np.e + 1.8 * beta * q) / (np.log(np.e + 1.8 * beta * q) + C * q**2)
        return T0
    f = 1.0 / (1.0 + (k * s / 5.4)**4)
    Tc = f * T_tilde(k, 1.0, beta_c) + (1.0 - f) * T_tilde(k, alpha_c, 1.0)
    y = z_eq / (1.0 + z_d)
    G = y * (-6.0 * np.sqrt(1.0 + y) + (2.0 + 3.0 * y) * np.log((np.sqrt(1.0 + y) + 1.0) / (np.sqrt(1.0 + y) - 1.0)))
    alpha_b = (2.07 * k_eq * s * (1.0 + R_d)**(-3.0 / 4.0) * G)
    beta_b = 0.5 + Omega_b / Omega_m + (3.0 - 2.0 * Omega_b / Omega_m) * np.sqrt((17.2 * ommh2)**2 + 1.0)
    beta_node = 8.41 * ommh2**0.435
    s_tilde = s / (1.0 + (beta_node / (k * s))**3)**(1.0 / 3.0)
    j0 = np.sinc(k * s_tilde / np.pi)
    Tb = (T_tilde(k, 1.0, 1.0) / (1.0 + (k * s / 5.2)**2) + alpha_b / (1.0 + (beta_b / (k * s))**3) * np.exp(-(k / k_silk)**1.4)) * j0
    T = Omega_b / Omega_m * Tb + (1.0 - Omega_b / Omega_m) * Tc
    return T

def tophat_window(x):
    return np.where(np.abs(x) < 1e-4, 1.0, 3.0 * (np.sin(x) - x * np.cos(x)) / x**3)

def virial_temperature_to_mass(T_vir_K, z):
    Hz = cosmo.H(z).to(u.km / u.s / u.Mpc).value * 1e3 / Mpc_m
    Delta_c = 18.0 * np.pi**2 + 82.0 * (cosmo.Om(z) - 1.0) - 39.0 * (cosmo.Om(z) - 1.0)**2
    M = (2.0 * kB * T_vir_K / (mu * mp))**(3.0 / 2.0) * (G_SI * Hz * np.sqrt(Delta_c) * (4.0 * np.pi / 3.0)**(1.0 / 3.0))**(-3.0 / 2.0)
    return M / Msun_kg

if __name__ == '__main__':
    z_vals = np.linspace(1, 20, 20)
    for z in [5, 10, 15, 20]:
        m_min = virial_temperature_to_mass(1e4, z)
        print('z=' + str(z) + ': M_min = ' + str(m_min) + ' M_sun')