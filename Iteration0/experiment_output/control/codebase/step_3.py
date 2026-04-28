# filename: codebase/step_3.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import stats
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
import os
import warnings

warnings.filterwarnings('ignore')

def load_results():
    theo = np.load('data/theoretical_model.npz')
    gems = np.load('data/gems_results.npy', allow_pickle=True).item()
    sparkler = np.load('data/sparkler_results.npy', allow_pickle=True).item()
    return theo, gems, sparkler

def build_theoretical_cdf(z_grid, rate, z_min, z_max):
    mask = (z_grid >= z_min) & (z_grid <= z_max)
    z_cdf = z_grid[mask]
    rate_cdf = rate[mask]
    rate_cdf_pos = np.clip(rate_cdf, 0.0, None)
    cum = cumulative_trapezoid(rate_cdf_pos, z_cdf, initial=0.0)
    total = cum[-1]
    if total <= 0.0:
        cdf_vals = np.linspace(0.0, 1.0, len(z_cdf))
    else:
        cdf_vals = cum / total
    cdf_interp = interp1d(z_cdf, cdf_vals, kind='linear', bounds_error=False, fill_value=(0.0, 1.0))
    return z_cdf, cdf_vals, cdf_interp

def ks_test_against_theoretical_cdf(median_zform, cdf_interp):
    sorted_z = np.sort(median_zform)
    n = len(sorted_z)
    ecdf_after = np.arange(1, n + 1) / n
    ecdf_before = np.arange(0, n) / n
    theo_cdf_vals = cdf_interp(sorted_z)
    d_plus = np.max(ecdf_after - theo_cdf_vals)
    d_minus = np.max(theo_cdf_vals - ecdf_before)
    ks_stat = max(d_plus, d_minus)
    p_value = stats.kstest(sorted_z, cdf_interp).pvalue
    return ks_stat, p_value

def peak_window_fraction(median_zform, z_grid, rate):
    peak_mask = rate >= 0.5
    z_in_window = z_grid[peak_mask]
    z_lo = z_in_window.min()
    z_hi = z_in_window.max()
    in_window = (median_zform >= z_lo) & (median_zform <= z_hi)
    n_in = int(np.sum(in_window))
    fraction = n_in / len(median_zform)
    return z_lo, z_hi, fraction, n_in

def partial_spearman_zform_zh_controlling_mass(zh, zform, log10_mass):
    n = len(zh)
    rank_zh = stats.rankdata(zh)
    rank_zform = stats.rankdata(zform)
    rank_mass = stats.rankdata(log10_mass)
    def ols_residuals(y_ranks, x_ranks):
        x_mean = np.mean(x_ranks)
        y_mean = np.mean(y_ranks)
        slope = np.sum((x_ranks - x_mean) * (y_ranks - y_mean)) / np.sum((x_ranks - x_mean) ** 2)
        intercept = y_mean - slope * x_mean
        fitted = slope * x_ranks + intercept
        return y_ranks - fitted
    resid_zh = ols_residuals(rank_zh, rank_mass)
    resid_zform = ols_residuals(rank_zform, rank_mass)
    partial_rho, _ = stats.pearsonr(resid_zh, resid_zform)
    df = n - 3
    if df <= 0:
        p_value = np.nan
    else:
        t_stat = partial_rho * np.sqrt(df) / np.sqrt(1.0 - partial_rho ** 2 + 1e-15)
        p_value = 2.0 * stats.t.sf(np.abs(t_stat), df=df)
    return partial_rho, p_value, n, df

if __name__ == '__main__':
    theo, gems, sparkler = load_results()
    z_grid = theo['z']
    rate = theo['rate']
    gems_median_zform = gems['zform_median']
    sparkler_median_zform = sparkler['zform_median']
    gems_zh = gems['ZH']
    sparkler_zh = sparkler['ZH']
    gems_mass = gems['log10_mass']
    sparkler_mass = sparkler['log10_mass']
    all_median_zform = np.concatenate([gems_median_zform, sparkler_median_zform])
    z_min_obs = float(np.min(all_median_zform))
    z_max_obs = float(np.max(all_median_zform))
    z_cdf, cdf_vals, cdf_interp = build_theoretical_cdf(z_grid, rate, z_min_obs, z_max_obs)
    ks_gems, p_gems = ks_test_against_theoretical_cdf(gems_median_zform, cdf_interp)
    ks_spark, p_spark = ks_test_against_theoretical_cdf(sparkler_median_zform, cdf_interp)
    print('GEMS KS: ' + str(ks_gems) + ', p: ' + str(p_gems))
    print('Sparkler KS: ' + str(ks_spark) + ', p: ' + str(p_spark))
    z_lo, z_hi, frac_g, _ = peak_window_fraction(gems_median_zform, z_grid, rate)
    z_lo, z_hi, frac_s, _ = peak_window_fraction(sparkler_median_zform, z_grid, rate)
    print('GEMS peak fraction: ' + str(frac_g))
    print('Sparkler peak fraction: ' + str(frac_s))
    rho_zh_zf, p_zh_zf = stats.spearmanr(np.concatenate([gems_zh, sparkler_zh]), np.concatenate([gems_median_zform, sparkler_median_zform]))
    print('Spearman [Z/H] vs z_form: rho=' + str(rho_zh_zf) + ', p=' + str(p_zh_zf))
    rho_zh_m, p_zh_m = stats.spearmanr(np.concatenate([gems_zh, sparkler_zh]), np.concatenate([gems_mass, sparkler_mass]))
    print('Spearman [Z/H] vs mass: rho=' + str(rho_zh_m) + ', p=' + str(p_zh_m))
    p_rho, p_pval, _, _ = partial_spearman_zform_zh_controlling_mass(np.concatenate([gems_zh, sparkler_zh]), np.concatenate([gems_median_zform, sparkler_median_zform]), np.concatenate([gems_mass, sparkler_mass]))
    print('Partial Spearman [Z/H] vs z_form (ctrl mass): rho=' + str(p_rho) + ', p=' + str(p_pval))