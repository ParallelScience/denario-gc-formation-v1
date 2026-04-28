# filename: codebase/step_3.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from scipy.stats import ks_2samp, anderson_ksamp, gaussian_kde
import warnings
import os

warnings.filterwarnings('ignore')

def compute_time_interval(cosmo, z1, z2):
    t_z1 = cosmo.age(z1).to(u.Gyr).value
    t_z2 = cosmo.age(z2).to(u.Gyr).value
    delta_t_gyr = t_z2 - t_z1
    return delta_t_gyr, t_z1, t_z2

def compute_survival_probabilities(log10_masses, alpha):
    masses = 10.0 ** log10_masses
    mass_max = np.max(masses)
    probs = (masses / mass_max) ** alpha
    return probs

def simulate_surviving_population(log10_masses, alpha, n_samples=200000, rng_seed=42):
    rng = np.random.default_rng(rng_seed)
    probs = compute_survival_probabilities(log10_masses, alpha)
    weights = probs / probs.sum()
    surviving_log10_masses = rng.choice(log10_masses, size=n_samples, replace=True, p=weights)
    return surviving_log10_masses, weights

def compute_kde_loglikelihood(observed_values, kde_samples, bandwidth_method='scott'):
    kde = gaussian_kde(kde_samples, bw_method=bandwidth_method)
    pdf_values = kde(observed_values)
    pdf_values = np.clip(pdf_values, 1e-300, None)
    log_likelihood = np.sum(np.log(pdf_values))
    likelihood = np.exp(log_likelihood)
    return log_likelihood, likelihood

def run_statistical_tests(simulated_masses, sparkler_masses):
    ks_stat, ks_pvalue = ks_2samp(simulated_masses, sparkler_masses)
    ad_result = anderson_ksamp([simulated_masses, sparkler_masses])
    ad_stat = ad_result.statistic
    ad_critical_values = ad_result.critical_values
    ad_significance_levels = ad_result.significance_level
    return ks_stat, ks_pvalue, ad_stat, ad_critical_values, ad_significance_levels

if __name__ == '__main__':
    data_dir = 'data/'
    cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315, Ob0=0.049)
    z_gems = 9.625
    z_sparkler = 1.4
    delta_t, t_gems, t_sparkler = compute_time_interval(cosmo, z_gems, z_sparkler)
    gc_df = pd.read_csv(os.path.join(data_dir, 'gc_data.csv'))
    gems_df = gc_df[gc_df['system'] == 'GEMS'].copy()
    sparkler_df = gc_df[gc_df['system'] == 'Sparkler'].copy()
    gems_masses = gems_df['log10_mass'].values
    sparkler_masses = sparkler_df['log10_mass'].values
    alpha_values = [0.5, 1.0, 2.0]
    n_samples = 200000
    results = []
    simulated_distributions = {}
    best_alpha = None
    best_loglik = -np.inf
    for alpha in alpha_values:
        surv_masses, weights = simulate_surviving_population(gems_masses, alpha, n_samples=n_samples, rng_seed=42)
        simulated_distributions[alpha] = surv_masses
        surv_mean = np.mean(surv_masses)
        surv_median = np.median(surv_masses)
        surv_std = np.std(surv_masses)
        ks_stat, ks_pvalue, ad_stat, ad_crit, ad_sig = run_statistical_tests(surv_masses, sparkler_masses)
        log_lik, lik = compute_kde_loglikelihood(sparkler_masses, surv_masses)
        if log_lik > best_loglik:
            best_loglik = log_lik
            best_alpha = alpha
        results.append({'alpha': alpha, 'surv_mean_log10M': surv_mean, 'surv_median_log10M': surv_median, 'surv_std_log10M': surv_std, 'ks_stat': ks_stat, 'ks_pvalue': ks_pvalue, 'ad_stat': ad_stat, 'log_likelihood': log_lik, 'likelihood': lik})
    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(data_dir, 'survival_modeling_results.csv'), index=False)
    save_dict = {}
    for alpha in alpha_values:
        save_dict['alpha_' + str(alpha).replace('.', 'p')] = simulated_distributions[alpha]
    save_dict['gems_log10_masses'] = gems_masses
    save_dict['sparkler_log10_masses'] = sparkler_masses
    save_dict['alpha_values'] = np.array(alpha_values)
    save_dict['best_alpha'] = np.array([best_alpha])
    save_dict['delta_t_gyr'] = np.array([delta_t])
    np.savez(os.path.join(data_dir, 'survival_simulated_distributions.npz'), **save_dict)