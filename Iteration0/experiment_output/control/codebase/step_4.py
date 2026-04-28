# filename: codebase/step_4.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, spearmanr, linregress
from astropy.cosmology import FlatLambdaCDM
import time
import warnings

warnings.filterwarnings('ignore')

def load_all_data():
    theo = np.load('data/theoretical_model.npz')
    gems = np.load('data/gems_results.npy', allow_pickle=True).item()
    sparkler = np.load('data/sparkler_results.npy', allow_pickle=True).item()
    return theo, gems, sparkler

def get_milky_way_gcs():
    mw_data = [(-9.19, -2.37), (-7.45, -2.27), (-7.37, -2.23), (-8.71, -2.10), (-8.93, -1.50), (-8.70, -1.53), (-8.20, -2.31), (-9.42, -0.72), (-9.41, -0.55), (-9.47, -0.44), (-7.49, -0.44), (-7.64, -0.64), (-7.73, -1.54), (-8.81, -1.29), (-7.32, -1.37)]
    mv_abs = np.array([d[0] for d in mw_data])
    mw_zh = np.array([d[1] for d in mw_data])
    mass_msun = 2.0 * 10.0 ** ((4.83 - mv_abs) / 2.5)
    return np.log10(mass_msun), mw_zh

if __name__ == '__main__':
    theo, gems, sparkler = load_all_data()
    cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315)
    ts = str(int(time.time()))
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(theo['z'], theo['rate'], color='black', lw=2, label='Atomic cooling model')
    ax.fill_between(theo['z'], 0, theo['rate'], where=theo['rate'] >= 0.5, alpha=0.25, color='gray', label='Peak window')
    ax.errorbar(gems['zform_median'], np.full(len(gems['zform_median']), 0.02), xerr=gems['zform_std'], fmt='o', color='royalblue', label='GEMS')
    ax.errorbar(sparkler['zform_median'], np.full(len(sparkler['zform_median']), 0.06), xerr=sparkler['zform_std'], fmt='s', color='tomato', label='Sparkler')
    ax.set_xlim(0, 32)
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('Normalized GC Formation Rate')
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    z_ticks = np.array([2, 5, 10, 15, 20, 25, 30])
    ax2.set_xticks(z_ticks)
    ax2.set_xticklabels([str(round(cosmo.lookback_time(z).value, 1)) for z in z_ticks], rotation=30, fontsize=8)
    ax2.set_xlabel('Lookback Time (Gyr)')
    ax.legend()
    plt.savefig('data/figure1_' + ts + '.png', dpi=300)
    plt.close()
    mw_mass, mw_zh = get_milky_way_gcs()
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(mw_mass, mw_zh, marker='^', color='green', label='MW GCs')
    ax.errorbar(gems['log10_mass'], gems['ZH'], yerr=[gems['ZH_err_lo'], gems['ZH_err_hi']], fmt='o', color='royalblue', label='GEMS')
    ax.errorbar(sparkler['log10_mass'], sparkler['ZH'], yerr=[sparkler['ZH_err_lo'], sparkler['ZH_err_hi']], fmt='s', color='tomato', label='Sparkler')
    all_m = np.concatenate([gems['log10_mass'], sparkler['log10_mass'], mw_mass])
    all_z = np.concatenate([gems['ZH'], sparkler['ZH'], mw_zh])
    slope, intercept, _, _, _ = linregress(all_m, all_z)
    x_line = np.linspace(all_m.min(), all_m.max(), 200)
    ax.plot(x_line, slope * x_line + intercept, color='black', linestyle='--', label='OLS Fit')
    ax.annotate('Slope: ' + str(round(slope, 2)) + ', Intercept: ' + str(round(intercept, 2)), xy=(0.05, 0.05), xycoords='axes fraction')
    ax.set_xlabel('log10(M/Msolar)')
    ax.set_ylabel('[Z/H]')
    ax.legend()
    plt.savefig('data/figure2_' + ts + '.png', dpi=300)
    plt.close()