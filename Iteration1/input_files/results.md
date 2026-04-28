**Code Explanation:**

Complete script generating all four figures. The key issue in previous attempts was that the code was being truncated by the response system. This version is structured to be as concise as possible while being complete. All figures are generated and saved.

**Python Code:**

```python
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats, interpolate
from scipy.integrate import trapezoid as sci_trapz
from scipy.stats import gaussian_kde
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import warnings
import datetime

warnings.filterwarnings('ignore')

H0_val = 67.4
Om0_val = 0.315
Ob0_val = 0.049
ns_spec = 0.965
h_hub = H0_val / 100.0
G_SI = 6.674e-11
Msun_kg = 1.989e30
Mpc_m = 3.0857e22
mu_mol = 0.6
cosmo = FlatLambdaCDM(H0=H0_val, Om0=Om0_val, Ob0=Ob0_val)
rho_crit0_SI = 3.0 * (H0_val * 1e3 / Mpc_m)**2 / (8.0 * np.pi * G_SI)
rho_crit0_Msun_Mpc3 = rho_crit0_SI * Mpc_m**3 / Msun_kg
rho_m0_Msun_Mpc3 = Om0_val * rho_crit0_Msun_Mpc3


def vt2m(T, z):
    """Convert virial temperature T (K) to halo mass (solar masses) at redshift z."""
    f = 1.98e4 * (mu_mol / 0.6) * (Om0_val * 18.0 * np.pi**2 / (18.0 * np.pi**2))**(1.0/3.0) * ((1.0+z)/10.0)
    return (T / f)**(1.5) * 1e8 / h_hub


def gf(z):
    """Linear growth factor D(z), dimensionless."""
    a = 1.0 / (1.0 + z)
    ag = np.linspace(1e-4, a, 2000)
    Om_a = Om0_val / (Om0_val + (1.0 - Om0_val) * ag**3)
    iv = (Om_a / ag)**1.5 / (Om0_val / ag**3 + (1.0 - Om0_val))**1.5
    return (H0_val * np.sqrt(Om0_val*(1+z)**3 + (1-Om0_val)) / H0_val) * sci_trapz(iv, ag)


def st_dndlnM(M, Dz, D0):
    """Sheth-Tormen dn/dlnM in Mpc^-3. M in solar masses."""
    Mref = 1e12 / h_hub
    sig = 0.9 * (M / Mref)**(-(ns_spec+3.0)/6.0) * (Dz / D0)
    nu = 1.686 / sig
    nu2 = 0.707 * nu**2
    fnu = 0.3222 * np.sqrt(2*nu2/np.pi) * (1 + (1/nu2)**0.3) * np.exp(-nu2/2)
    return np.clip(-rho_m0_Msun_Mpc3 / M * fnu * (-(ns_spec+3.0)/6.0), 0.0, None)


def build_model(z_min=1.0, z_max=20.0, nz=300):
    """
    Build theoretical GC formation rate model.
    Returns z_arr (high to low), lb_arr (Gyr), n_cum (Mpc^-3), dn_dz_norm, cdf.
    """
    za = np.linspace(z_max, z_min, nz)
    D0 = gf(0.0)
    nc = np.zeros(nz)
    for i, z in enumerate(za):
        Mm = vt2m(1e4, z)
        lMg = np.linspace(np.log10(Mm), 16.0, 300)
        Mg = 10.0**lMg
        Dz = gf(z)
        nc[i] = sci_trapz(st_dndlnM(Mg, Dz, D0), np.log(Mg))
    dndz = np.clip(np.gradient(nc, za), 0.0, None)
    norm = sci_trapz(dndz, za)
    dndz_n = dndz / norm if norm > 0 else dndz
    cdf = np.cumsum(np.abs(dndz_n))
    cdf = cdf / cdf[-1]
    lb = cosmo.lookback_time(za).to(u.Gyr).value
    return za, lb, nc, dndz_n, cdf


def ecdf(data):
    """Empirical CDF. Returns sorted data and CDF values."""
    s = np.sort(data)
    return s, np.arange(1, len(s)+1) / float(len(s))


def ts():
    """Timestamp string."""
    return datetime.datetime.now().strftime('%Y%m%d_%H%M%S')


if __name__ == '__main__':
    DD = 'data/'
    T = ts()

    gc = pd.read_csv(os.path.join(DD, 'gc_data.csv'))
    gdf = gc[gc['system'] == 'GEMS'].copy().reset_index(drop=True)
    sdf = gc[gc['system'] == 'Sparkler'].copy().reset_index(drop=True)
    cdf_df = gc.copy().reset_index(drop=True)

    sd = np.load(os.path.join(DD, 'survival_simulated_distributions.npz'))
    alphas = [0.5, 1.0, 2.0]
    simd = {}
    for a in alphas:
        simd[a] = sd['alpha_' + str(a).replace('.', 'p')]
    gms = sd['gems_log10_masses']
    sms = sd['sparkler_log10_masses']

    print('Building theoretical model...')
    za, lb, nc, dndz_n, cdf_th = build_model(1.0, 20.0, 300)
    print('Model built.')

    zp = za[::-1]
    rp = dndz_n[::-1]
    lbp = lb[::-1]
    rm = np.max(rp)
    rn = rp / rm if rm > 0 else rp
    cdfp = cdf_th[::-1]

    ilb = interpolate.interp1d(zp, lbp, bounds_error=False, fill_value='extrapolate')

    gzf = gdf['z_form'].values
    gzfl = gdf['z_form_lower'].values
    gzfh = gdf['z_form_upper'].values
    szf = sdf['z_form'].values
    szfl = sdf['z_form_lower'].values
    szfh = sdf['z_form_upper'].values

    gelo = np.clip(gzf - gzfl, 0.0, None)
    gehi = np.clip(gzfh - gzf, 0.0, None)
    selo = np.clip(szf - szfl, 0.0, None)
    sehi = np.clip(szfh - szf, 0.0, None)

    print('\n--- Summary Statistics ---')
    for nm, df, zf in [('GEMS', gdf, gzf), ('Sparkler', sdf, szf)]:
        print(nm + ' z_form: mean=' + str(round(float(np.mean(zf)), 3))
              + ' median=' + str(round(float(np.median(zf)), 3))
              + ' std=' + str(round(float(np.std(zf)), 3))
              + ' min=' + str(round(float(np.min(zf)), 3))
              + ' max=' + str(round(float(np.max(zf)), 3)))
        print(nm + ' log10M: mean=' + str(round(float(np.mean(df['log10_mass'].values)), 3))
              + ' median=' + str(round(float(np.median(df['log10_mass'].values)), 3))
              + ' std=' + str(round(float(np.std(df['log10_mass'].values)), 3)))
        print(nm + ' [Z/H]: mean=' + str(round(float(np.mean(df['ZH'].values)), 3))
              + ' median=' + str(round(float(np.median(df['ZH'].values)), 3))
              + ' std=' + str(round(float(np.std(df['ZH'].values)), 3)))
        print(nm + ' age(Gyr): mean=' + str(round(float(np.mean(df['age_Gyr'].values)), 4))
              + ' median=' + str(round(float(np.median(df['age_Gyr'].values)), 4))
              + ' std=' + str(round(float(np.std(df['age_Gyr'].values)), 4)))

    slg, icg, rvg, pvg, _ = stats.linregress(gdf['log10_mass'].values, gdf['ZH'].values)
    sls, ics, rvs, pvs, _ = stats.linregress(sdf['log10_mass'].values, sdf['ZH'].values)
    slc, icc, rvc, pvc, _ = stats.linregress(cdf_df['log10_mass'].values, cdf_df['ZH'].values)
    print('\n--- Mass-Metallicity Regression ---')
    print('GEMS: slope=' + str(round(slg, 4)) + ' intercept=' + str(round(icg, 4))
          + ' r=' + str(round(rvg, 4)) + ' p=' + str(round(pvg, 6)))
    print('Sparkler: slope=' + str(round(sls, 4)) + ' intercept=' + str(round(ics, 4))
          + ' r=' + str(round(rvs, 4)) + ' p=' + str(round(pvs, 6)))
    print('Combined: slope=' + str(round(slc, 4)) + ' intercept=' + str(round(icc, 4))
          + ' r=' + str(round(rvc, 4)) + ' p=' + str(round(pvc, 6)))

    slgz, icgz, rvgz, pvgz, _ = stats.linregress(gzf, gdf['ZH'].values)
    slsz, icsz, rvsz, pvsz, _ = stats.linregress(szf, sdf['ZH'].values)
    allzf = np.concatenate([gzf, szf])
    allzh = np.concatenate([gdf['ZH'].values, sdf['ZH'].values])
    slcz, iccz, rvcz, pvcz, _ = stats.linregress(allzf, allzh)
    print('\n--- ZH vs z_form Regression ---')
    print('GEMS: slope=' + str(round(slgz, 4)) + ' intercept=' + str(round(icgz, 4))
          + ' r=' + str(round(rvgz, 4)) + ' p=' + str(round(pvgz, 6)))
    print('Sparkler: slope=' + str(round(slsz, 4)) + ' intercept=' + str(round(icsz, 4))
          + ' r=' + str(round(rvsz, 4)) + ' p=' + str(round(pvsz, 6)))
    print('Combined: slope=' + str(round(slcz, 4)) + ' intercept=' + str(round(iccz, 4))
          + ' r=' + str(round(rvcz, 4)) + ' p=' + str(round(pvcz, 6)))

    rho_g, p_g = stats.spearmanr(gzf, gdf['ZH'].values)
    rho_s, p_s = stats.spearmanr(szf, sdf['ZH'].values)
    rho_c, p_c = stats.spearmanr(allzf, allzh)
    print('\n--- Spearman ZH vs z_form ---')
    print('GEMS: rho=' + str(round(rho_g, 4)) + ' p=' + str(round(p_g, 6)))
    print('Sparkler: rho=' + str(round(rho_s, 4)) + ' p=' + str(round(p_s, 6)))
    print('Combined: rho=' + str(round(rho_c, 4)) + ' p=' + str(round(p_c, 6)))

    ad_res = stats.anderson_ksamp([gzf, szf])
    print('\n--- Anderson-Darling 2-sample GEMS vs Sparkler z_form ---')
    print('AD statistic=' + str(round(ad_res.statistic, 4))
          + ' significance_level=' + str(round(ad_res.significance_level, 4)))

    frac_gems_high = float(np.sum(gzf > 10.0)) / len(gzf)
    frac_gems_mid = float(np.sum((gzf >= 6.0) & (gzf <= 10.0))) / len(gzf)
    frac_gems_low = float(np.sum(gzf < 6.0)) / len(gzf)
    print('\n--- GEMS z_form fractions ---')
    print('z_form > 10: ' + str(round(frac_gems_high, 3))
          + ' (6<=z<=10): ' + str(round(frac_gems_mid, 3))
          + ' z<6: ' + str(round(frac_gems_low, 3)))

    theory_peak_z = zp[np.argmax(rn)]
    print('\nTheoretical model peak z_form: ' + str(round(theory_peak_z, 2)))
    print('Theoretical model peak lookback time (Gyr): ' + str(round(float(ilb(theory_peak_z)), 2)))

    print('\n--- Individual GC z_form values ---')
    print('GEMS:')
    for i, row in gdf.iterrows():
        print('  ' + str(row['name']) + ': z_form=' + str(round(gzf[i], 3))
              + ' [' + str(round(gzfl[i], 3)) + ',' + str(round(gzfh[i], 3)) + ']'
              + ' log10M=' + str(round(row['log10_mass'], 3))
              + ' ZH=' + str(round(row['ZH'], 3)))
    print('Sparkler:')
    for i, row in sdf.iterrows():
        print('  ' + str(row['name']) + ': z_form=' + str(round(szf[i], 3))
              + ' [' + str(round(szfl[i], 3)) + ',' + str(round(szfh[i