# filename: codebase/step_3.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
from scipy import stats, interpolate

def load_data(data_dir):
    df = pd.read_csv(os.path.join(data_dir, 'gc_formation_redshifts.csv'))
    npz = np.load(os.path.join(data_dir, 'trenti_model.npz'))
    model = {'z': npz['z'], 'pdf': npz['pdf'], 'cdf': npz['cdf']}
    return df, model

def build_theoretical_cdf_func(model):
    z_arr = model['z']
    cdf_arr = model['cdf']
    cdf_func = interpolate.interp1d(z_arr, cdf_arr, kind='linear', bounds_error=False, fill_value=(0.0, 1.0))
    return cdf_func

def ks_one_sample_gems(gems_zform, cdf_func):
    ks_stat, p_value = stats.kstest(gems_zform, cdf_func)
    return ks_stat, p_value

def ks_two_sample(gems_zform, sparkler_zform):
    ks_stat, p_value = stats.ks_2samp(gems_zform, sparkler_zform)
    return ks_stat, p_value

def spearman_zh_zform(zh_arr, zform_arr, label):
    rho, p_value = stats.spearmanr(zh_arr, zform_arr)
    return rho, p_value

def compute_summary_table(df):
    cols = ['age_Gyr', 'log10_mass', 'ZH', 'z_form']
    rows = []
    for system in ['GEMS', 'Sparkler', 'Combined']:
        if system == 'Combined':
            sub = df[cols].dropna()
        else:
            sub = df[df['system'] == system][cols].dropna()
        for stat_name, func in [('mean', np.mean), ('median', np.median), ('std', np.std)]:
            row = {'system': system, 'statistic': stat_name}
            for c in cols:
                row[c] = func(sub[c].values)
            rows.append(row)
    summary = pd.DataFrame(rows)
    summary = summary.set_index(['system', 'statistic'])
    return summary

if __name__ == '__main__':
    pd.set_option('display.max_rows', 200)
    pd.set_option('display.max_columns', 20)
    pd.set_option('display.width', 200)
    pd.set_option('display.float_format', '{:.6f}'.format)
    data_dir = 'data/'
    df, model = load_data(data_dir)
    gems_df = df[df['system'] == 'GEMS'].dropna(subset=['z_form'])
    sparkler_df = df[df['system'] == 'Sparkler'].dropna(subset=['z_form'])
    gems_zform = gems_df['z_form'].values
    sparkler_zform = sparkler_df['z_form'].values
    cdf_func = build_theoretical_cdf_func(model)
    ks1_stat, ks1_pval = ks_one_sample_gems(gems_zform, cdf_func)
    ks2_stat, ks2_pval = ks_two_sample(gems_zform, sparkler_zform)
    combined_df = df.dropna(subset=['z_form'])
    rho_comb, p_comb = spearman_zh_zform(combined_df['ZH'].values, combined_df['z_form'].values, 'Combined')
    rho_gems, p_gems = spearman_zh_zform(gems_df['ZH'].values, gems_zform, 'GEMS')
    rho_spark, p_spark = spearman_zh_zform(sparkler_df['ZH'].values, sparkler_zform, 'Sparkler')
    summary = compute_summary_table(df)
    print('=' * 70)
    print('ONE-SAMPLE K-S TEST: GEMS z_form vs. Trenti et al. (2015) CDF')
    print('  NOTE: Limited sample size N=' + str(len(gems_zform)) + ' reduces statistical power.')
    print('  K-S statistic : ' + str(round(ks1_stat, 6)))
    print('  p-value        : ' + str(round(ks1_pval, 6)))
    print('=' * 70)
    print('TWO-SAMPLE K-S TEST: GEMS z_form vs. Sparkler z_form')
    print('  NOTE: Very small Sparkler sample N=' + str(len(sparkler_zform)) + '; results are indicative only.')
    print('  K-S statistic : ' + str(round(ks2_stat, 6)))
    print('  p-value        : ' + str(round(ks2_pval, 6)))
    print('=' * 70)
    print('SPEARMAN RANK CORRELATION: [Z/H] vs. z_form')
    print('  Combined (N=' + str(len(combined_df)) + '): rho = ' + str(round(rho_comb, 6)) + ', p-value = ' + str(round(p_comb, 6)))
    print('  GEMS only (N=' + str(len(gems_df)) + '): rho = ' + str(round(rho_gems, 6)) + ', p-value = ' + str(round(p_gems, 6)))
    print('  Sparkler only (N=' + str(len(sparkler_df)) + '): rho = ' + str(round(rho_spark, 6)) + ', p-value = ' + str(round(p_spark, 6)))
    print('=' * 70)
    print('STATISTICAL SUMMARY TABLE')
    print(summary.to_string())
    stat_results = pd.DataFrame([{'test': 'KS_one_sample_GEMS_vs_Trenti', 'statistic': ks1_stat, 'p_value': ks1_pval, 'note': 'N=19 GEMS z_form vs theoretical CDF'}, {'test': 'KS_two_sample_GEMS_vs_Sparkler', 'statistic': ks2_stat, 'p_value': ks2_pval, 'note': 'N=19 GEMS vs N=5 Sparkler z_form'}, {'test': 'Spearman_ZH_zform_Combined', 'statistic': rho_comb, 'p_value': p_comb, 'note': 'N=' + str(len(combined_df)) + ' combined sample'}, {'test': 'Spearman_ZH_zform_GEMS', 'statistic': rho_gems, 'p_value': p_gems, 'note': 'N=' + str(len(gems_df)) + ' GEMS only'}, {'test': 'Spearman_ZH_zform_Sparkler', 'statistic': rho_spark, 'p_value': p_spark, 'note': 'N=' + str(len(sparkler_df)) + ' Sparkler only'}])
    stat_results.to_csv(os.path.join(data_dir, 'statistical_results.csv'), index=False)
    summary.to_csv(os.path.join(data_dir, 'summary_statistics.csv'))
    print('Statistical results saved to data/statistical_results.csv')
    print('Summary statistics saved to data/summary_statistics.csv')