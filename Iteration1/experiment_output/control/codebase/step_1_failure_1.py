# filename: codebase/step_1.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.cosmology import z_at_value
import os
import warnings
warnings.filterwarnings('ignore')

def compute_z_form(cosmo, z_obs, age_gyr, age_min_gyr=0.01):
    t_obs = cosmo.age(z_obs).to(u.Gyr).value
    t_form = t_obs - age_gyr
    t_form_clamped = max(t_form, age_min_gyr)
    z_form = z_at_value(cosmo.age, t_form_clamped * u.Gyr, zmin=0.01, zmax=50.0).value
    return z_form, t_form_clamped

def build_gc_dataframe():
    gems_data = [
        ("J1", 7.679883, 0.306606, 0.105830, 0.156109, -0.155737, 0.172521, 0.425764),
        ("I1", 7.209972, 0.222194, 0.118820, 0.103790, -0.116777, 0.168080, 0.353989),
        ("H1", 7.112587, 0.268571, 0.158832, 0.136833, -0.166758, 0.206958, 0.414991),
        ("G1", 6.580490, 0.091426, 0.103001, 0.064732, -0.178940, 0.213542, 0.486345),
        ("A1", 7.489207, 0.108196, 0.047334, 0.040569, -0.204631, 0.218476, 0.424188),
        ("B1", 8.047918, 0.207793, 0.072150, 0.071099, -0.005409, 0.082002, 0.163142),
        ("C1", 7.517093, 0.111555, 0.067088, 0.063398, -0.060425, 0.128039, 0.266884),
        ("D1", 7.882942, 0.220086, 0.097852, 0.100520, -0.127497, 0.175482, 0.291695),
        ("E1", 6.766869, 0.056349, 0.075598, 0.039125, -0.184665, 0.220689, 0.572239),
        ("F1", 6.909087, 0.097596, 0.083108, 0.058072, -0.159828, 0.200231, 0.463759),
        ("F2", 7.091763, 0.132655, 0.056626, 0.061231, -0.271403, 0.272512, 0.729343),
        ("E2", 6.202470, 0.003594, 0.002773, 0.001775, -0.179677, 0.221968, 0.346192),
        ("D2", 7.327680, 0.113312, 0.093868, 0.076996, -0.166798, 0.205681, 0.485723),
        ("C2", 7.575265, 0.115583, 0.093113, 0.073473, -0.109043, 0.168116, 0.337416),
        ("B2", 7.915161, 0.149981, 0.047094, 0.045620, 0.041284, 0.051739, 0.115378),
        ("A2", 7.223565, 0.090363, 0.053208, 0.034850, -0.609474, 0.464776, 0.440869),
        ("H2", 6.643740, 0.091809, 0.118365, 0.058557, -0.236564, 0.258077, 0.635067),
        ("I2", 7.273337, 0.247759, 0.168230, 0.142788, -0.115895, 0.152703, 0.495809),
        ("J2", 7.155033, 0.192443, 0.088012, 0.097319, -0.154601, 0.201931, 0.464464),
    ]
    gems_df = pd.DataFrame(gems_data, columns=["name", "log10_mass", "age_Gyr", "age_err_lo", "age_err_hi", "ZH", "ZH_err_lo", "ZH_err_hi"])
    gems_df["system"] = "GEMS"
    gems_df["z_obs"] = 9.625
    sparkler_data = [
        ("SP1", 6.891181, 2.730595, 5.221295, 1.528180, -0.849845, 0.658266, 0.529515),
        ("SP2", 7.026021, 2.028381, 1.500056, 1.044356, -0.482711, 0.347040, 0.568824),
        ("SP4", 6.846949, 1.508617, 1.239569, 0.714472, -0.403906, 0.315182, 0.608286),
        ("SP8", 7.065139, 1.743875, 1.094816, 0.785341, -0.964864, 0.542399, 0.566541),
        ("SP10", 6.944087, 1.660122, 1.029555, 0.969398, -0.185263, 0.190706, 0.366281),
    ]
    sparkler_df = pd.DataFrame(sparkler_data, columns=["name", "log10_mass", "age_Gyr", "age_err_lo", "age_err_hi", "ZH", "ZH_err_lo", "ZH_err_hi"])
    sparkler_df["system"] = "Sparkler"
    sparkler_df["z_obs"] = 1.4
    return gems_df, sparkler_df

def process_gc_dataframe(df, cosmo):
    z_form_list, z_form_lower_list, z_form_upper_list, t_form_list, lookback_list = [], [], [], [], []
    for _, row in df.iterrows():
        z_obs, age, age_lo, age_hi = row["z_obs"], row["age_Gyr"], row["age_err_lo"], row["age_err_hi"]
        zf, tf = compute_z_form(cosmo, z_obs, age)
        zf_lower, _ = compute_z_form(cosmo, z_obs, age + age_hi)
        zf_upper, _ = compute_z_form(cosmo, z_obs, age - age_lo)
        lb = cosmo.lookback_time(zf).to(u.Gyr).value
        z_form_list.append(zf); z_form_lower_list.append(zf_lower); z_form_upper_list.append(zf_upper); t_form_list.append(tf); lookback_list.append(lb)
    df = df.copy()
    df["z_form"] = z_form_list; df["z_form_lower"] = z_form_lower_list; df["z_form_upper"] = z_form_upper_list; df["t_form_gyr"] = t_form_list; df["lookback_time_form_Gyr"] = lookback_list
    return df

if __name__ == '__main__':
    cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315)
    gems_df, sparkler_df = build_gc_dataframe()
    gems_df = process_gc_dataframe(gems_df, cosmo)
    sparkler_df = process_gc_dataframe(sparkler_df, cosmo)
    all_df = pd.concat([gems_df, sparkler_df])
    all_df.to_csv("data/gc_data.csv", index=False)
    print("Saved to data/gc_data.csv")