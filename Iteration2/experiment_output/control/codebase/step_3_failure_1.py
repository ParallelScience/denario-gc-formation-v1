# filename: codebase/step_3.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats, integrate
from scipy.stats import theilslopes, spearmanr
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import os
import time
from step_1 import get_cosmology, tvir_to_mhalo

def sheth_tormen_f(nu):
    A = 0.3222
    a = 0.707
    p = 0.3
    nu2 = nu ** 2
    return A * np.sqrt(2.0 * a * nu2 / np.pi) * (1.0 + (a * nu2) ** (-p)) * np.exp(-0.5 * a * nu2)

def growth_factor(z, cosmo):
    def integrand(zp):
        Ez = cosmo.efunc(zp)
        return (1.0 + zp) / Ez ** 3
    def D_unnorm(zz):
        val, _ = integrate.quad(integrand, zz, 1000.0)
        return cosmo.efunc(zz) * val
    return D_unnorm(z) / D_unnorm(0.0)

def sigma_m(M_msun, z, cosmo, sigma8=0.811, ns=0.965):
    M8 = (4.0 / 3.0) * np.pi * (8.0 ** 3) * cosmo.critical_density(0).to(u.Msun / u.Mpc**3).value * cosmo.Om0
    n_eff = ns - 1.0
    sigma0 = sigma8 * (M_msun / M8) ** (-(3.0 + n_eff + 3.0) / 6.0)
    Dz = growth_factor(z, cosmo)
    return sigma0 * Dz

def dn_dlnM(M_msun, z, cosmo, sigma8=0.811, ns=0.965, delta_c=1.686):
    rho_m = (cosmo.Om0 * cosmo.critical_density(0).to(u.Msun / u.Mpc**3).value)
    sig = sigma_m(M_msun, z, cosmo, sigma8, ns)
    nu = delta_c / sig
    dln_sigma_dlnM = -(3.0 + ns - 1.0 + 3.0) / 6.0
    f_nu = sheth_tormen_f(nu)
    return (rho_m / M_msun) * f_nu * nu * abs(dln_sigma_dlnM)

def compute_n_above_mmin(z, cosmo, T_vir=1e4, sigma8=0.811, ns=0.965, n_pts=60):
    M_min = tvir_to_mhalo(T_vir, z, cosmo)
    M_max = 1e15
    if M_min >= M_max: return 0.0
    lnM_arr = np.linspace(np.log(M_min), np.log(M_max), n_pts)
    M_arr = np.exp(lnM_arr)
    integrand = np.array([dn_dlnM(m, z, cosmo, sigma8, ns) for m in M_arr])
    return np.trapz(integrand, lnM_arr)

def compute_dn_dz(z_arr, cosmo, sigma8=0.811, ns=0.965):
    n_arr = np.array([compute_n_above_mmin(z, cosmo, sigma8=sigma8, ns=ns) for z in z_arr])
    dn_dz = np.abs(np.gradient(n_arr, z_arr))
    return dn_dz, n_arr

if __name__ == '__main__':
    data_dir = 'data/'
    cosmo = get_cosmology()
    df = pd.read_csv(os.path.join(data_dir, 'zform_all_gcs.csv'))
    df_gems = df[df['system'] == 'GEMS'].copy()
    df_sparkler = df[df['system'] == 'Sparkler'].copy()
    z_theory = np.linspace(0.5, 20.0, 80)
    dn_dz_vals, n_vals = compute_dn_dz(z_theory, cosmo)
    dn_dz_norm = dn_dz_vals / np.nanmax(dn_dz_vals)
    fig = plt.figure(figsize=(18, 14))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.35)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.plot(z_theory, dn_dz_norm, 'k-', lw=2)
    for _, row in df_gems.iterrows():
        ax_a.axvline(row['z_form'], color='#1f77b4', alpha=0.35, lw=0.8)
    for _, row in df_sparkler.iterrows():
        ax_a.axvline(row['z_form'], color='#ff7f0e', alpha=0.5, lw=1.0)
    plt.savefig(os.path.join(data_dir, 'analysis_plot_' + str(int(time.time())) + '.png'))
    print('Saved to ' + os.path.join(data_dir, 'analysis_plot_' + str(int(time.time())) + '.png'))