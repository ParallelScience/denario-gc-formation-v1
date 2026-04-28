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
from scipy.integrate import trapezoid
from scipy.stats import theilslopes, spearmanr
import astropy.units as u
import time
from step_1 import get_cosmology, tvir_to_mhalo

def sheth_tormen_f(nu):
    A, a, p = 0.3222, 0.707, 0.3
    nu2 = nu ** 2
    return A * np.sqrt(2.0 * a * nu2 / np.pi) * (1.0 + (a * nu2) ** (-p)) * np.exp(-0.5 * a * nu2)

def growth_factor(z, cosmo):
    def integrand(zp):
        return (1.0 + zp) / cosmo.efunc(zp) ** 3
    def D_unnorm(zz):
        val, _ = integrate.quad(integrand, zz, 1000.0)
        return cosmo.efunc(zz) * val
    return D_unnorm(z) / D_unnorm(0.0)

def compute_n_above_mmin(z, cosmo, rho_m, n_eff, sigma8=0.811, ns=0.965, n_pts=60, delta_c=1.686, T_vir=1e4):
    M_min = tvir_to_mhalo(T_vir, z, cosmo)
    M_max = 1e15
    if M_min >= M_max: return 0.0
    M8 = (4.0 / 3.0) * np.pi * (8.0 ** 3) * rho_m
    Dz = growth_factor(z, cosmo)
    lnM_arr = np.linspace(np.log(M_min), np.log(M_max), n_pts)
    M_arr = np.exp(lnM_arr)
    dln_sigma_dlnM = -(3.0 + n_eff + 3.0) / 6.0
    sigma0_arr = sigma8 * (M_arr / M8) ** dln_sigma_dlnM
    sig_arr = sigma0_arr * Dz
    nu_arr = delta_c / sig_arr
    f_nu_arr = sheth_tormen_f(nu_arr)
    integrand_arr = (rho_m / M_arr) * f_nu_arr * nu_arr * abs(dln_sigma_dlnM)
    return trapezoid(integrand_arr, lnM_arr)

def compute_dn_dz_curve(z_arr, cosmo, sigma8=0.811, ns=0.965):
    rho_m = cosmo.Om0 * cosmo.critical_density(0).to(u.Msun / u.Mpc**3).value
    n_eff = ns - 1.0
    n_arr = np.array([compute_n_above_mmin(z, cosmo, rho_m, n_eff, sigma8=sigma8, ns=ns) for z in z_arr])
    dn_dz = np.abs(np.gradient(n_arr, z_arr))
    return dn_dz, n_arr

def ccmf_fit(log10_masses):
    masses = 10.0 ** np.array(log10_masses)
    sorted_m = np.sort(masses)
    n = len(sorted_m)
    ccdf = np.arange(n, 0, -1, dtype=float)
    slope, intercept, r, p, se = stats.linregress(np.log10(sorted_m), np.log10(ccdf))
    return sorted_m, ccdf, slope, intercept, 1.0 - slope

if __name__ == '__main__':
    data_dir = 'data/'
    cosmo = get_cosmology()
    df = pd.read_csv(os.path.join(data_dir, 'zform_all_gcs.csv'))
    df_gems = df[df['system'] == 'GEMS'].copy().reset_index(drop=True)
    df_sparkler = df[df['system'] == 'Sparkler'].copy().reset_index(drop=True)
    z_theory = np.linspace(0.5, 20.0, 80)
    dn_dz_vals, n_vals = compute_dn_dz_curve(z_theory, cosmo)
    dn_dz_norm = dn_dz_vals / np.nanmax(dn_dz_vals)
    fig, axes = plt.subplots(3, 2, figsize=(16, 18))
    ax_a, ax_b, ax_c, ax_d, ax_e, ax_f = axes.flatten()
    for ax, lbl in zip(axes.flatten(), ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']): ax.text(0.02, 0.97, lbl, transform=ax.transAxes, fontsize=13, fontweight='bold')
    ax_a.plot(z_theory, dn_dz_norm, 'k-', lw=2.0)
    ax_b.plot(z_theory, dn_dz_norm, 'k-', lw=2.0)
    ax_c.scatter(df_gems['log10_mass'], df_gems['ZH'], color='#1f77b4', label='GEMS')
    ax_c.scatter(df_sparkler['log10_mass'], df_sparkler['ZH'], color='#ff7f0e', label='Sparkler')
    ax_d.scatter(df_gems['z_form'], df_gems['ZH'], color='#1f77b4')
    ax_d.scatter(df_sparkler['z_form'], df_sparkler['ZH'], color='#ff7f0e')
    ax_e.loglog(np.sort(10**df_gems['log10_mass']), np.arange(len(df_gems), 0, -1), 'o', color='#1f77b4')
    ax_f.hist(df_gems['age_Gyr'], alpha=0.5, color='#1f77b4', label='GEMS')
    ax_f.hist(df_sparkler['age_Gyr'], alpha=0.5, color='#ff7f0e', label='Sparkler')
    plt.tight_layout()
    plt.savefig(os.path.join(data_dir, 'analysis_plot_' + str(int(time.time())) + '.png'))
    print('Saved to ' + os.path.join(data_dir, 'analysis_plot_' + str(int(time.time())) + '.png'))