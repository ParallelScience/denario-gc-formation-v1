# filename: codebase/step_2.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
import os
import warnings

warnings.filterwarnings('ignore')

H0_val = 67.4
Om0_val = 0.315
Ob0_val = 0.049
cosmo = FlatLambdaCDM(H0=H0_val, Om0=Om0_val, Ob0=Ob0_val)

SENTINEL_Z = 50.0
N_MC = 1000
RNG_SEED = 42

GEMS_DATA = [
    ('J1', 0.0, 7.679883, 0.306606, 0.105830, 0.156109, -0.155737, 0.172521, 0.425764),
    ('I1', 0.0, 7.209972, 0.222194, 0.118820, 0.103790, -0.116777, 0.168080, 0.353989),
    ('H1', 0.0, 7.112587, 0.268571, 0.158832, 0.136833, -0.166758, 0.206958, 0.414991),
    ('G1', 0.0, 6.580490, 0.091426, 0.103001, 0.064732, -0.178940, 0.213542, 0.486345),
    ('A1', 0.0, 7.489207, 0.108196, 0.047334, 0.040569, -0.204631, 0.218476, 0.424188),
    ('B1', 0.0, 8.047918, 0.207793, 0.072150, 0.071099, -0.005409, 0.082002, 0.163142),
    ('C1', 0.0, 7.517093, 0.111555, 0.067088, 0.063398, -0.060425, 0.128039, 0.266884),
    ('D1', 0.0, 7.882942, 0.220086, 0.097852, 0.100520, -0.127497, 0.175482, 0.291695),
    ('E1', 0.0, 6.766869, 0.056349, 0.075598, 0.039125, -0.184665, 0.220689, 0.572239),
    ('F1', 0.0, 6.909087, 0.097596, 0.083108, 0.058072, -0.159828, 0.200231, 0.463759),
    ('F2', 0.0, 7.091763, 0.132655, 0.056626, 0.061231, -0.271403, 0.272512, 0.729343),
    ('E2', 0.0, 6.202470, 0.003594, 0.002773, 0.001775, -0.179677, 0.221968, 0.346192),
    ('D2', 0.0, 7.327680, 0.113312, 0.093868, 0.076996, -0.166798, 0.205681, 0.485723),
    ('C2', 0.0, 7.575265, 0.115583, 0.093113, 0.073473, -0.109043, 0.168116, 0.337416),
    ('B2', 0.0, 7.915161, 0.149981, 0.047094, 0.045620, 0.041284, 0.051739, 0.115378),
    ('A2', 0.0, 7.223565, 0.090363, 0.053208, 0.034850, -0.609474, 0.464776, 0.440869),
    ('H2', 0.0, 6.643740, 0.091809, 0.118365, 0.058557, -0.236564, 0.258077, 0.635067),
    ('I2', 0.0, 7.273337, 0.247759, 0.168230, 0.142788, -0.115895, 0.152703, 0.495809),
    ('J2', 0.0, 7.155033, 0.192443, 0.088012, 0.097319, -0.154601, 0.201931, 0.464464),
]

SPARKLER_DATA = [
    (1, 6.891181, 2.730595, 5.221295, 1.528180, -0.849845, 0.658266, 0.529515),
    (2, 7.026021, 2.028381, 1.500056, 1.044356, -0.482711, 0.347040, 0.568824),
    (4, 6.846949, 1.508617, 1.239569, 0.714472, -0.403906, 0.315182, 0.608286),
    (8, 7.065139, 1.743875, 1.094816, 0.785341, -0.964864, 0.542399, 0.566541),
    (10, 6.944087, 1.660122, 1.029555, 0.969398, -0.185263, 0.190706, 0.366281),
]

Z_OBS_GEMS = 9.625
Z_OBS_SPARKLER = 1.4

def build_z_from_lookback_interpolator(cosmo_obj, z_max=60.0, n_points=10000):
    z_arr = np.linspace(0.0, z_max, n_points)
    lt_arr = cosmo_obj.lookback_time(z_arr).value
    return interp1d(lt_arr, z_arr, kind='linear', bounds_error=False, fill_value=(0.0, z_max))

def age_universe_at_z(z, cosmo_obj):
    return cosmo_obj.age(z).value

def compute_zform_point(age_gyr, z_obs, cosmo_obj, lt_to_z_interp):
    lt_obs = cosmo_obj.lookback_time(z_obs).value
    lt_form = lt_obs + age_gyr
    return float(lt_to_z_interp(lt_form))

def mc_zform_samples(age_gyr, age_err_lo, age_err_hi, z_obs, cosmo_obj, lt_to_z_interp, n_samples, rng):
    age_max = age_universe_at_z(z_obs, cosmo_obj)
    samples = np.empty(n_samples)
    filled = 0
    batch = max(n_samples * 4, 4000)
    while filled < n_samples:
        u_vals = rng.uniform(0.0, 1.0, batch)
        gauss = rng.standard_normal(batch)
        sigma = np.where(gauss < 0.0, age_err_lo, age_err_hi)
        age_draws = age_gyr + sigma * gauss
        valid = (age_draws > 0.0) & (age_draws < age_max)
        valid_draws = age_draws[valid]
        n_take = min(len(valid_draws), n_samples - filled)
        samples[filled:filled + n_take] = valid_draws[:n_take]
        filled += n_take
    lt_obs = cosmo_obj.lookback_time(z_obs).value
    lt_form_samples = lt_obs + samples
    return lt_to_z_interp(lt_form_samples).astype(float)

def process_system(data_list, z_obs, system_name, cosmo_obj, lt_to_z_interp, rng):
    age_universe_obs = age_universe_at_z(z_obs, cosmo_obj)
    names, log10_mass_list, age_list, ZH_list, ZH_err_lo_list, ZH_err_hi_list, zform_point_list, zform_mc_list = [], [], [], [], [], [], [], []
    warned = False
    for row in data_list:
        if system_name == 'GEMS':
            name, ebv, lmass, age, err_lo, err_hi, zh, zh_lo, zh_hi = row
        else:
            name, lmass, age, err_lo, err_hi, zh, zh_lo, zh_hi = row
        if age >= age_universe_obs:
            if not warned:
                print('WARNING: ' + system_name + ' GC ' + str(name) + ' age=' + str(round(age, 4)) + ' Gyr >= age of universe at z_obs=' + str(round(age_universe_obs, 4)) + ' Gyr. Capping z_form at sentinel=' + str(SENTINEL_Z))
                warned = True
            zform_pt = SENTINEL_Z
        else:
            zform_pt = compute_zform_point(age, z_obs, cosmo_obj, lt_to_z_interp)
        zform_mc = mc_zform_samples(age, err_lo, err_hi, z_obs, cosmo_obj, lt_to_z_interp, N_MC, rng)
        names.append(str(name))
        log10_mass_list.append(lmass)
        age_list.append(age)
        ZH_list.append(zh)
        ZH_err_lo_list.append(zh_lo)
        ZH_err_hi_list.append(zh_hi)
        zform_point_list.append(zform_pt)
        zform_mc_list.append(zform_mc)
    zform_mc_arr = np.array(zform_mc_list)
    return {'names': names, 'log10_mass': np.array(log10_mass_list), 'age_gyr': np.array(age_list), 'ZH': np.array(ZH_list), 'ZH_err_lo': np.array(ZH_err_lo_list), 'ZH_err_hi': np.array(ZH_err_hi_list), 'zform_point': np.array(zform_point_list), 'zform_mc': zform_mc_arr, 'zform_mean': np.mean(zform_mc_arr, axis=1), 'zform_median': np.median(zform_mc_arr, axis=1), 'zform_std': np.std(zform_mc_arr, axis=1)}

def print_per_gc_stats(result, system_name):
    print('\n--- ' + system_name + ' Per-GC Statistics ---')
    for i, name in enumerate(result['names']):
        print('GC ' + name + ': z_form_mean=' + str(round(result['zform_mean'][i], 3)) + ', z_form_std=' + str(round(result['zform_std'][i], 3)))

if __name__ == '__main__':
    rng = np.random.default_rng(RNG_SEED)
    lt_to_z = build_z_from_lookback_interpolator(cosmo)
    gems_res = process_system(GEMS_DATA, Z_OBS_GEMS, 'GEMS', cosmo, lt_to_z, rng)
    sparkler_res = process_system(SPARKLER_DATA, Z_OBS_SPARKLER, 'Sparkler', cosmo, lt_to_z, rng)
    print_per_gc_stats(gems_res, 'GEMS')
    print_per_gc_stats(sparkler_res, 'Sparkler')
    np.save('data/gems_results.npy', gems_res)
    np.save('data/sparkler_results.npy', sparkler_res)
    print('\nResults saved to data/ directory.')