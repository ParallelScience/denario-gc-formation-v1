# filename: codebase/step_5.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
import os
from step_1 import get_cosmology, tvir_to_mhalo

def stellar_to_halo_mass(log10_mstar, epsilon):
    f_b = 0.049 / 0.315
    log10_mhalo = log10_mstar - np.log10(epsilon * f_b)
    return log10_mhalo

def maiolino2008_mzr(log10_mstar, z):
    z_nodes = np.array([0.07, 0.7, 2.2, 3.5])
    logM0_nodes = np.array([11.18, 11.57, 12.38, 12.76])
    K0_nodes = np.array([9.185, 9.185, 8.971, 8.901])
    z_clipped = np.clip(z, z_nodes[0], z_nodes[-1])
    logM0 = np.interp(z_clipped, z_nodes, logM0_nodes)
    K0 = np.interp(z_clipped, z_nodes, K0_nodes)
    oh = -0.0864 * (log10_mstar - logM0) ** 2 + K0
    return oh

def oh_to_zh(oh, oh_solar=8.69):
    return oh - oh_solar

if __name__ == '__main__':
    data_dir = 'data/'
    cosmo = get_cosmology()
    df = pd.read_csv(os.path.join(data_dir, 'zform_all_gcs.csv'))
    df_gems = df[df['system'] == 'GEMS'].copy().reset_index(drop=True)
    df_sparkler = df[df['system'] == 'Sparkler'].copy().reset_index(drop=True)
    f_b = 0.049 / 0.315
    epsilon_lo = 0.01
    epsilon_hi = 0.05
    print('=' * 70)
    print('STEP 5: MOLECULAR COOLING THRESHOLD AND MZR ANALYSIS')
    print('=' * 70)
    print('\n--- Stellar-to-Halo Mass Conversion (Sparkler GCs) ---')
    print('Assumed epsilon range: ' + str(epsilon_lo) + ' to ' + str(epsilon_hi))
    print('Baryon fraction f_b = Omega_b/Omega_m = 0.049/0.315 = ' + '%.4f' % f_b)
    print('Justification: At high redshift (z~2-4), low-mass halos (M_halo~10^7-10^8 Msolar) have shallow potential wells and are subject to strong supernova feedback and photoionization heating, suppressing star formation efficiency to epsilon~1-5%.')
    print('\nSparkler GC implied halo masses:')
    print('  GC_ID  log10(M_star)  z_form  log10(M_halo,eps=0.01)  log10(M_halo,eps=0.05)')
    sparkler_halo_lo = []
    sparkler_halo_hi = []
    for _, row in df_sparkler.iterrows():
        mh_lo = stellar_to_halo_mass(row['log10_mass'], epsilon_hi)
        mh_hi = stellar_to_halo_mass(row['log10_mass'], epsilon_lo)
        sparkler_halo_lo.append(mh_lo)
        sparkler_halo_hi.append(mh_hi)
        print('  ' + str(row['name']) + '  ' + '%.4f' % row['log10_mass'] + '  ' + '%.3f' % row['z_form'] + '  ' + '%.3f' % mh_hi + '  ' + '%.3f' % mh_lo)
    sparkler_halo_lo = np.array(sparkler_halo_lo)
    sparkler_halo_hi = np.array(sparkler_halo_hi)
    print('\nSparkler implied halo mass range:')
    print('  epsilon=0.05: log10(M_halo) = ' + '%.3f' % np.min(sparkler_halo_lo) + ' to ' + '%.3f' % np.max(sparkler_halo_lo) + ' (Msolar)')
    print('  epsilon=0.01: log10(M_halo) = ' + '%.3f' % np.min(sparkler_halo_hi) + ' to ' + '%.3f' % np.max(sparkler_halo_hi) + ' (Msolar)')
    mol_cool_lo = 5.0
    mol_cool_hi = 6.0
    atom_cool_lo = 7.0
    atom_cool_hi = 8.0
    print('\n--- Cooling Threshold Comparison ---')
    print('Molecular cooling threshold: log10(M_halo) ~ ' + str(mol_cool_lo) + ' to ' + str(mol_cool_hi) + ' (Msolar)')
    print('Atomic cooling threshold:    log10(M_halo) ~ ' + str(atom_cool_lo) + ' to ' + str(atom_cool_hi) + ' (Msolar)')
    median_halo_lo = np.median(sparkler_halo_lo)
    median_halo_hi = np.median(sparkler_halo_hi)
    in_atomic_lo = np.sum((sparkler_halo_lo >= atom_cool_lo) & (sparkler_halo_lo <= atom_cool_hi))
    in_atomic_hi = np.sum((sparkler_halo_hi >= atom_cool_lo) & (sparkler_halo_hi <= atom_cool_hi))
    above_atomic = np.sum(sparkler_halo_lo > atom_cool_hi)
    print('\nNumber of Sparkler GCs with implied halo mass in atomic cooling range [7,8]:')
    print('  epsilon=0.05: ' + str(in_atomic_lo) + ' / ' + str(len(df_sparkler)))
    print('  epsilon=0.01: ' + str(in_atomic_hi) + ' / ' + str(len(df_sparkler)))
    if median_halo_lo >= atom_cool_lo:
        verdict = 'CONSISTENT WITH ATOMIC COOLING HALOS (T_vir >= 10^4 K)'
    elif median_halo_hi >= atom_cool_lo:
        verdict = 'MARGINALLY CONSISTENT WITH ATOMIC COOLING HALOS'
    else:
        verdict = 'CONSISTENT WITH MOLECULAR COOLING HALOS (T_vir < 10^4 K)'
    print('\nVerdict: Sparkler GCs are ' + verdict)
    print('\n--- Maiolino et al. (2008) MZR Analysis ---')
    median_logm_gems = np.median(df_gems['log10_mass'].values)
    median_zform_gems = np.median(df_gems['z_form'].dropna().values)
    median_zh_gems_obs = np.median(df_gems['ZH'].values)
    median_logm_sparkler = np.median(df_sparkler['log10_mass'].values)
    median_zform_sparkler = np.median(df_sparkler['z_form'].dropna().values)
    median_zh_sparkler_obs = np.median(df_sparkler['ZH'].values)
    oh_gems_pred = maiolino2008_mzr(median_logm_gems, median_zform_gems)
    zh_gems_pred = oh_to_zh(oh_gems_pred)
    oh_sparkler_pred = maiolino2008_mzr(median_logm_sparkler, median_zform_sparkler)
    zh_sparkler_pred = oh_to_zh(oh_sparkler_pred)
    print('\nGEMS system:')
    print('  Median log10(M_star/Msolar) = ' + '%.4f' % median_logm_gems)
    print('  Median z_form               = ' + '%.4f' % median_zform_gems)
    print('  MZR predicted [Z/H]         = ' + '%.4f' % zh_gems_pred)
    print('  Observed median [Z/H]       = ' + '%.4f' % median_zh_gems_obs)
    print('\nSparkler system:')
    print('  Median log10(M_star/Msolar) = ' + '%.4f' % median_logm_sparkler)
    print('  Median z_form               = ' + '%.4f' % median_zform_sparkler)
    print('  MZR predicted [Z/H]         = ' + '%.4f' % zh_sparkler_pred)
    print('  Observed median [Z/H]       = ' + '%.4f' % median_zh_sparkler_obs)