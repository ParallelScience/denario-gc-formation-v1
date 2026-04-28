# filename: codebase/step_1.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.cosmology import z_at_value

def compute_z_form(lookback_form_Gyr, cosmo, z_min=0.001, z_max=50.0):
    age_universe = cosmo.age(0).to(u.Gyr).value
    if lookback_form_Gyr <= 0.0 or lookback_form_Gyr >= age_universe:
        return np.nan
    try:
        z_f = z_at_value(cosmo.lookback_time, lookback_form_Gyr * u.Gyr, zmin=z_min, zmax=z_max)
        return float(z_f)
    except Exception:
        return np.nan

def build_gc_table(cosmo):
    gems_names = ['J1','I1','H1','G1','A1','B1','C1','D1','E1','F1','F2','E2','D2','C2','B2','A2','H2','I2','J2']
    gems_data = np.array([
        [7.679883, 0.306606, 0.105830, 0.156109, -0.155737, 0.172521, 0.425764],
        [7.209972, 0.222194, 0.118820, 0.103790, -0.116777, 0.168080, 0.353989],
        [7.112587, 0.268571, 0.158832, 0.136833, -0.166758, 0.206958, 0.414991],
        [6.580490, 0.091426, 0.103001, 0.064732, -0.178940, 0.213542, 0.486345],
        [7.489207, 0.108196, 0.047334, 0.040569, -0.204631, 0.218476, 0.424188],
        [8.047918, 0.207793, 0.072150, 0.071099, -0.005409, 0.082002, 0.163142],
        [7.517093, 0.111555, 0.067088, 0.063398, -0.060425, 0.128039, 0.266884],
        [7.882942, 0.220086, 0.097852, 0.100520, -0.127497, 0.175482, 0.291695],
        [6.766869, 0.056349, 0.075598, 0.039125, -0.184665, 0.220689, 0.572239],
        [6.909087, 0.097596, 0.083108, 0.058072, -0.159828, 0.200231, 0.463759],
        [7.091763, 0.132655, 0.056626, 0.061231, -0.271403, 0.272512, 0.729343],
        [6.202470, 0.003594, 0.002773, 0.001775, -0.179677, 0.221968, 0.346192],
        [7.327680, 0.113312, 0.093868, 0.076996, -0.166798, 0.205681, 0.485723],
        [7.575265, 0.115583, 0.093113, 0.073473, -0.109043, 0.168116, 0.337416],
        [7.915161, 0.149981, 0.047094, 0.045620,  0.041284, 0.051739, 0.115378],
        [7.223565, 0.090363, 0.053208, 0.034850, -0.609474, 0.464776, 0.440869],
        [6.643740, 0.091809, 0.118365, 0.058557, -0.236564, 0.258077, 0.635067],
        [7.273337, 0.247759, 0.168230, 0.142788, -0.115895, 0.152703, 0.495809],
        [7.155033, 0.192443, 0.088012, 0.097319, -0.154601, 0.201931, 0.464464]
    ])
    sparkler_names = ['S1','S2','S4','S8','S10']
    sparkler_data = np.array([
        [6.891181, 2.730595, 5.221295, 1.528180, -0.849845, 0.658266, 0.529515],
        [7.026021, 2.028381, 1.500056, 1.044356, -0.482711, 0.347040, 0.568824],
        [6.846949, 1.508617, 1.239569, 0.714472, -0.403906, 0.315182, 0.608286],
        [7.065139, 1.743875, 1.094816, 0.785341, -0.964864, 0.542399, 0.566541],
        [6.944087, 1.660122, 1.029555, 0.969398, -0.185263, 0.190706, 0.366281]
    ])
    z_gems, z_sparkler = 9.625, 1.4
    lt_gems, lt_sparkler = cosmo.lookback_time(z_gems).to(u.Gyr).value, cosmo.lookback_time(z_sparkler).to(u.Gyr).value
    rows = []
    for i, name in enumerate(gems_names):
        log10_mass, age, age_err_lo, age_err_hi, ZH, ZH_err_lo, ZH_err_hi = gems_data[i]
        lt_form_med, lt_form_lo, lt_form_hi = lt_gems - age, lt_gems - (age + age_err_hi), lt_gems - (age - age_err_lo)
        rows.append({'system': 'GEMS', 'name': name, 'z_obs': z_gems, 'age_Gyr': age, 'age_err_lo': age_err_lo, 'age_err_hi': age_err_hi, 'log10_mass': log10_mass, 'ZH': ZH, 'ZH_err_lo': ZH_err_lo, 'ZH_err_hi': ZH_err_hi, 'z_form': compute_z_form(lt_form_med, cosmo), 'z_form_lower_bound': compute_z_form(lt_form_lo, cosmo), 'z_form_upper_bound': compute_z_form(lt_form_hi, cosmo)})
    for i, name in enumerate(sparkler_names):
        log10_mass, age, age_err_lo, age_err_hi, ZH, ZH_err_lo, ZH_err_hi = sparkler_data[i]
        lt_form_med, lt_form_lo, lt_form_hi = lt_sparkler - age, lt_sparkler - (age + age_err_hi), lt_sparkler - (age - age_err_lo)
        rows.append({'system': 'Sparkler', 'name': name, 'z_obs': z_sparkler, 'age_Gyr': age, 'age_err_lo': age_err_lo, 'age_err_hi': age_err_hi, 'log10_mass': log10_mass, 'ZH': ZH, 'ZH_err_lo': ZH_err_lo, 'ZH_err_hi': ZH_err_hi, 'z_form': compute_z_form(lt_form_med, cosmo), 'z_form_lower_bound': compute_z_form(lt_form_lo, cosmo), 'z_form_upper_bound': compute_z_form(lt_form_hi, cosmo)})
    return pd.DataFrame(rows)

if __name__ == '__main__':
    cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315, Ob0=0.049)
    df = build_gc_table(cosmo)
    df.to_csv('data/gc_formation_redshifts.csv', index=False)
    print('Saved to data/gc_formation_redshifts.csv')