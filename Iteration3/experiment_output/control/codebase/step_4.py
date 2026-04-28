# filename: codebase/step_4.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.cosmology import z_at_value
import datetime

def compute_z_form_correct(z_obs, age_Gyr, cosmo, z_max=50.0):
    age_universe = cosmo.age(0).to(u.Gyr).value
    lt_obs = cosmo.lookback_time(z_obs).to(u.Gyr).value
    lt_form = lt_obs + age_Gyr
    if lt_form <= 0.0 or lt_form >= age_universe:
        return np.nan
    try:
        z_f = z_at_value(cosmo.lookback_time, lt_form * u.Gyr, zmin=z_obs, zmax=z_max)
        return float(z_f)
    except Exception:
        return np.nan

if __name__ == '__main__':
    data_dir = 'data/'
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315, Ob0=0.049)
    df = pd.read_csv(os.path.join(data_dir, 'gc_formation_redshifts.csv'))
    npz = np.load(os.path.join(data_dir, 'trenti_model.npz'))
    stat_results = pd.read_csv(os.path.join(data_dir, 'statistical_results.csv'))
    z_model, pdf_model, cdf_model = npz['z'], npz['pdf'], npz['cdf']
    sparkler_log10_mass = np.array([6.891181, 7.026021, 6.846949, 7.065139, 6.944087])
    sparkler_age_Gyr = np.array([2.730595, 2.028381, 1.508617, 1.743875, 1.660122])
    sparkler_age_err_lo = np.array([5.221295, 1.500056, 1.239569, 1.094816, 1.029555])
    sparkler_age_err_hi = np.array([1.528180, 1.044356, 0.714472, 0.785341, 0.969398])
    sparkler_ZH = np.array([-0.849845, -0.482711, -0.403906, -0.964864, -0.185263])
    sparkler_ZH_err_lo = np.array([0.658266, 0.347040, 0.315182, 0.542399, 0.190706])
    sparkler_ZH_err_hi = np.array([0.529515, 0.568824, 0.608286, 0.566541, 0.366281])
    z_sparkler_obs = 1.4
    spark_zform = np.array([compute_z_form_correct(z_sparkler_obs, a, cosmo) for a in sparkler_age_Gyr])
    spark_zform_lo = np.array([compute_z_form_correct(z_sparkler_obs, a + e, cosmo) for a, e in zip(sparkler_age_Gyr, sparkler_age_err_hi)])
    spark_zform_hi = np.array([compute_z_form_correct(z_sparkler_obs, max(a - e, 0.001), cosmo) for a, e in zip(sparkler_age_Gyr, sparkler_age_err_lo)])
    gems_df = df[df['system'] == 'GEMS'].dropna(subset=['z_form']).copy()
    gems_zform, gems_zform_lo_v, gems_zform_hi_v = gems_df['z_form'].values, gems_df['z_form_lower_bound'].values, gems_df['z_form_upper_bound'].values
    gems_zh, gems_zh_lo, gems_zh_hi, gems_mass = gems_df['ZH'].values, gems_df['ZH_err_lo'].values, gems_df['ZH_err_hi'].values, gems_df['log10_mass'].values
    gems_xerr_lo, gems_xerr_hi = np.clip(np.where(np.isnan(gems_zform_lo_v), 0.0, gems_zform - gems_zform_lo_v), 0.0, None), np.clip(np.where(np.isnan(gems_zform_hi_v), 0.0, gems_zform_hi_v - gems_zform), 0.0, None)
    valid_s = ~np.isnan(spark_zform)
    spark_xerr_lo, spark_xerr_hi = np.clip(np.where(np.isnan(spark_zform_lo), 0.0, spark_zform - spark_zform_lo), 0.0, None), np.clip(np.where(np.isnan(spark_zform_hi), 0.0, spark_zform_hi - spark_zform), 0.0, None)
    ks2_pval = float(stat_results[stat_results['test'] == 'KS_two_sample_GEMS_vs_Sparkler'].iloc[0]['p_value'])
    gems_color, sparkler_color, model_color = '#1f77b4', '#d62728', '#2ca02c'
    fig1, axes1 = plt.subplots(2, 1, figsize=(11, 10))
    axes1[0].fill_between(z_model, pdf_model, alpha=0.35, color=model_color, label='Trenti+2015 model PDF')
    axes1[0].errorbar(gems_zform, np.full(len(gems_zform), np.max(pdf_model)*0.85), xerr=[gems_xerr_lo, gems_xerr_hi], fmt='o', color=gems_color, label='GEMS')
    axes1[0].errorbar(spark_zform[valid_s], np.full(valid_s.sum(), np.max(pdf_model)*0.7), xerr=[spark_xerr_lo[valid_s], spark_xerr_hi[valid_s]], fmt='s', color=sparkler_color, label='Sparkler')
    axes1[0].legend(); axes1[0].grid(True, alpha=0.3)
    axes1[1].plot(z_model, cdf_model, color=model_color, lw=2.5, label='Trenti+2015 model CDF')
    axes1[1].step(np.sort(gems_zform), np.arange(1, len(gems_zform)+1)/len(gems_zform), where='post', color=gems_color, label='GEMS')
    axes1[1].step(np.sort(spark_zform[valid_s]), np.arange(1, valid_s.sum()+1)/valid_s.sum(), where='post', color=sparkler_color, ls='--', label='Sparkler')
    axes1[1].legend(); axes1[1].grid(True, alpha=0.3)
    fig1.savefig(os.path.join(data_dir, 'formation_redshift_comparison_1_' + timestamp + '.png'), dpi=300)
    fig2, axes2 = plt.subplots(1, 3, figsize=(18, 6))
    axes2[0].errorbar(gems_mass, gems_zh, yerr=[gems_zh_lo, gems_zh_hi], fmt='o', color=gems_color, label='GEMS')
    axes2[0].errorbar(sparkler_log10_mass, sparkler_ZH, yerr=[sparkler_ZH_err_lo, sparkler_ZH_err_hi], fmt='s', color=sparkler_color, label='Sparkler')
    axes2[0].plot(np.linspace(5.5, 9.0, 100), 0.4*np.linspace(5.5, 9.0, 100)-3.0, 'k--', label='Harris+2016 MZR')
    axes2[0].legend(); axes2[0].grid(True, alpha=0.3)
    axes2[1].errorbar(gems_zform, gems_zh, xerr=[gems_xerr_lo, gems_xerr_hi], yerr=[gems_zh_lo, gems_zh_hi], fmt='o', color=gems_color, label='GEMS')
    axes2[1].errorbar(spark_zform[valid_s], sparkler_ZH[valid_s], xerr=[spark_xerr_lo[valid_s], spark_xerr_hi[valid_s]], yerr=[sparkler_ZH_err_lo[valid_s], sparkler_ZH_err_hi[valid_s]], fmt='s', color=sparkler_color, label='Sparkler')
    axes2[1].legend(); axes2[1].grid(True, alpha=0.3)
    for data, color, label in [(gems_mass, gems_color, 'GEMS'), (sparkler_log10_mass, sparkler_color, 'Sparkler')]:
        kde = gaussian_kde(data, bw_method=0.4)
        axes2[2].plot(np.linspace(5.5, 9.0, 300), kde(np.linspace(5.5, 9.0, 300)), color=color, label=label)
    axes2[2].text(0.03, 0.97, 'p={:.2e}'.format(ks2_pval), transform=axes2[2].transAxes, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    axes2[2].legend(); axes2[2].grid(True, alpha=0.3)
    fig2.savefig(os.path.join(data_dir, 'combined_analysis_2_' + timestamp + '.png'), dpi=300)