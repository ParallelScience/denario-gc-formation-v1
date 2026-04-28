# filename: codebase/step_2.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import theilslopes, ks_2samp, spearmanr
import os
from step_1 import get_cosmology, compute_zform_with_errors
if __name__ == '__main__':
    data_dir = 'data/'
    cosmo = get_cosmology()
    gems_data = [('J1', 0.0, 7.679883, 0.306606, 0.105830, 0.156109, -0.155737, 0.172521, 0.425764), ('I1', 0.0, 7.209972, 0.222194, 0.118820, 0.103790, -0.116777, 0.168080, 0.353989), ('H1', 0.0, 7.112587, 0.268571, 0.158832, 0.136833, -0.166758, 0.206958, 0.414991), ('G1', 0.0, 6.580490, 0.091426, 0.103001, 0.064732, -0.178940, 0.213542, 0.486345), ('A1', 0.0, 7.489207, 0.108196, 0.047334, 0.040569, -0.204631, 0.218476, 0.424188), ('B1', 0.0, 8.047918, 0.207793, 0.072150, 0.071099, -0.005409, 0.082002, 0.163142), ('C1', 0.0, 7.517093, 0.111555, 0.067088, 0.063398, -0.060425, 0.128039, 0.266884), ('D1', 0.0, 7.882942, 0.220086, 0.097852, 0.100520, -0.127497, 0.175482, 0.291695), ('E1', 0.0, 6.766869, 0.056349, 0.075598, 0.039125, -0.184665, 0.220689, 0.572239), ('F1', 0.0, 6.909087, 0.097596, 0.083108, 0.058072, -0.159828, 0.200231, 0.463759), ('F2', 0.0, 7.091763, 0.132655, 0.056626, 0.061231, -0.271403, 0.272512, 0.729343), ('E2', 0.0, 6.202470, 0.003594, 0.002773, 0.001775, -0.179677, 0.221968, 0.346192), ('D2', 0.0, 7.327680, 0.113312, 0.093868, 0.076996, -0.166798, 0.205681, 0.485723), ('C2', 0.0, 7.575265, 0.115583, 0.093113, 0.073473, -0.109043, 0.168116, 0.337416), ('B2', 0.0, 7.915161, 0.149981, 0.047094, 0.045620, 0.041284, 0.051739, 0.115378), ('A2', 0.0, 7.223565, 0.090363, 0.053208, 0.034850, -0.609474, 0.464776, 0.440869), ('H2', 0.0, 6.643740, 0.091809, 0.118365, 0.058557, -0.236564, 0.258077, 0.635067), ('I2', 0.0, 7.273337, 0.247759, 0.168230, 0.142788, -0.115895, 0.152703, 0.495809), ('J2', 0.0, 7.155033, 0.192443, 0.088012, 0.097319, -0.154601, 0.201931, 0.464464)]
    sparkler_data = [(1, 6.891181, 2.730595, 5.221295, 1.528180, -0.849845, 0.658266, 0.529515), (2, 7.026021, 2.028381, 1.500056, 1.044356, -0.482711, 0.347040, 0.568824), (4, 6.846949, 1.508617, 1.239569, 0.714472, -0.403906, 0.315182, 0.608286), (8, 7.065139, 1.743875, 1.094816, 0.785341, -0.964864, 0.542399, 0.566541), (10, 6.944087, 1.660122, 1.029555, 0.969398, -0.185263, 0.190706, 0.366281)]
    z_obs_gems, z_obs_sparkler = 9.625, 1.4
    gems_rows = []
    for row in gems_data:
        name, ebv, logm, age, age_lo, age_hi, zh, zh_lo, zh_hi = row
        zc, ze_lo, ze_hi = compute_zform_with_errors(z_obs_gems, age, age_lo, age_hi, cosmo)
        gems_rows.append({'name': name, 'system': 'GEMS', 'z_obs': z_obs_gems, 'log10_mass': logm, 'age_Gyr': age, 'age_err_lo': age_lo, 'age_err_hi': age_hi, 'ZH': zh, 'ZH_err_lo': zh_lo, 'ZH_err_hi': zh_hi, 'z_form': zc, 'z_form_err_lo': ze_lo, 'z_form_err_hi': ze_hi})
    sparkler_rows = []
    for row in sparkler_data:
        gid, logm, age, age_lo, age_hi, zh, zh_lo, zh_hi = row
        zc, ze_lo, ze_hi = compute_zform_with_errors(z_obs_sparkler, age, age_lo, age_hi, cosmo)
        sparkler_rows.append({'name': str(gid), 'system': 'Sparkler', 'z_obs': z_obs_sparkler, 'log10_mass': logm, 'age_Gyr': age, 'age_err_lo': age_lo, 'age_err_hi': age_hi, 'ZH': zh, 'ZH_err_lo': zh_lo, 'ZH_err_hi': zh_hi, 'z_form': zc, 'z_form_err_lo': ze_lo, 'z_form_err_hi': ze_hi})
    df_gems, df_sparkler = pd.DataFrame(gems_rows), pd.DataFrame(sparkler_rows)
    df_all = pd.concat([df_gems, df_sparkler], ignore_index=True)
    df_all.to_csv(os.path.join(data_dir, 'zform_all_gcs.csv'), index=False)
    def print_stats(df, label):
        print('\n--- Statistical Summary: ' + label + ' (N=' + str(len(df)) + ') ---')
        for col, unit in [('age_Gyr', 'Gyr'), ('log10_mass', 'log10(M/Msolar)'), ('ZH', '[Z/H]'), ('z_form', 'dimensionless')]:
            vals = df[col].dropna().values
            print('  ' + col + ' [' + unit + ']: mean=' + '%.4f' % np.mean(vals) + '  median=' + '%.4f' % np.median(vals) + '  std=' + '%.4f' % np.std(vals, ddof=1) + '  min=' + '%.4f' % np.min(vals) + '  max=' + '%.4f' % np.max(vals))
    print_stats(df_gems, 'GEMS')
    print_stats(df_sparkler, 'Sparkler')
    def theilsen_regression(x, y, label):
        res = theilslopes(y, x, 0.95)
        print('\n  Theil-Sen [Z/H] vs log10(M) [' + label + ']:')
        print('    slope=' + '%.4f' % res.slope + '  intercept=' + '%.4f' % res.intercept + '  slope_95CI=[' + '%.4f' % res.low_slope + ', ' + '%.4f' % res.high_slope + ']')
        return res.slope, res.intercept
    theilsen_regression(df_gems['log10_mass'].values, df_gems['ZH'].values, 'GEMS')
    theilsen_regression(df_sparkler['log10_mass'].values, df_sparkler['ZH'].values, 'Sparkler')
    theilsen_regression(df_all['log10_mass'].values, df_all['ZH'].values, 'Combined')
    def mass_function_slope(log10_masses, label):
        masses = 10.0 ** np.array(log10_masses)
        sorted_m = np.sort(masses)
        n = len(sorted_m)
        ccdf = np.arange(n, 0, -1, dtype=float)
        slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(sorted_m), np.log10(ccdf))
        print('\n  Mass function slope [' + label + ']: alpha (dN/dM)=' + '%.4f' % (1.0 - slope))
    mass_function_slope(df_gems['log10_mass'].values, 'GEMS')
    mass_function_slope(df_sparkler['log10_mass'].values, 'Sparkler')
    ks_stat, ks_pval = ks_2samp(df_gems['log10_mass'].values, df_sparkler['log10_mass'].values)
    print('\n=== K-S Test: GEMS vs Sparkler Mass Distributions ===')
    print('  KS stat=' + '%.4f' % ks_stat + ', p-value=' + '%.4f' % ks_pval)
    rho, p = spearmanr(df_all['ZH'], df_all['z_form'])
    print('\n=== Spearman Correlation: [Z/H] vs z_form (Combined) ===')
    print('  rho=' + '%.4f' % rho + ', p-value=' + '%.4f' % p)