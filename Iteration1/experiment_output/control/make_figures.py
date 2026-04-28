import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = False
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats, interpolate
from scipy.integrate import trapezoid as sci_trapz
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import warnings
warnings.filterwarnings('ignore')

H0_val = 67.4; Om0_val = 0.315; Ob0_val = 0.049; ns_spec = 0.965
h_hub = H0_val / 100.0
G_SI = 6.674e-11; Msun_kg = 1.989e30; Mpc_m = 3.0857e22; mu_mol = 0.6
cosmo = FlatLambdaCDM(H0=H0_val, Om0=Om0_val, Ob0=Ob0_val)
rho_crit0_SI = 3.0*(H0_val*1e3/Mpc_m)**2/(8.0*np.pi*G_SI)
rho_crit0 = rho_crit0_SI * Mpc_m**3 / Msun_kg
rho_m0 = Om0_val * rho_crit0

def vt2m(T, z):
    f = 1.98e4*(mu_mol/0.6)*(Om0_val*18.0*np.pi**2/(18.0*np.pi**2))**(1.0/3.0)*((1.0+z)/10.0)
    return (T/f)**1.5 * 1e8 / h_hub

def gf(z):
    a = 1.0/(1.0+z)
    ag = np.linspace(1e-4, a, 2000)
    Om_a = Om0_val/(Om0_val+(1.0-Om0_val)*ag**3)
    iv = (Om_a/ag)**1.5/(Om0_val/ag**3+(1.0-Om0_val))**1.5
    return (H0_val*np.sqrt(Om0_val*(1+z)**3+(1-Om0_val))/H0_val)*sci_trapz(iv, ag)

def st_dndlnM(M, Dz, D0):
    Mref = 1e12/h_hub
    sig = 0.9*(M/Mref)**(-(ns_spec+3.0)/6.0)*(Dz/D0)
    nu = 1.686/sig
    nu2 = 0.707*nu**2
    fnu = 0.3222*np.sqrt(2*nu2/np.pi)*(1+(1/nu2)**0.3)*np.exp(-nu2/2)
    return np.clip(-rho_m0/M*fnu*(-(ns_spec+3.0)/6.0), 0.0, None)

def build_model(z_min=1.0, z_max=20.0, nz=300):
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
    dndz_n = dndz/norm if norm > 0 else dndz
    cdf = np.cumsum(np.abs(dndz_n)); cdf = cdf/cdf[-1]
    lb = cosmo.lookback_time(za).to(u.Gyr).value
    return za, lb, nc, dndz_n, cdf

print('Building theoretical model...')
za, lb, nc, dndz_n, cdf_th = build_model(1.0, 20.0, 300)
zp = za[::-1]; rp = dndz_n[::-1]; lbp = lb[::-1]
nc_p = nc[::-1]
rm = np.max(rp); rn = rp/rm if rm > 0 else rp
cdfp = cdf_th[::-1]
ilb = interpolate.interp1d(zp, lbp, bounds_error=False, fill_value='extrapolate')
print('Model built.')

DD = 'data/'
gc = pd.read_csv(DD+'gc_data.csv')
gdf = gc[gc['system']=='GEMS'].copy().reset_index(drop=True)
sdf = gc[gc['system']=='Sparkler'].copy().reset_index(drop=True)
surv = pd.read_csv(DD+'survival_modeling_results.csv')
sd = np.load(DD+'survival_simulated_distributions.npz')

gzf = gdf['z_form'].values; gzfl = gdf['z_form_lower'].values; gzfh = gdf['z_form_upper'].values
szf = sdf['z_form'].values; szfl = sdf['z_form_lower'].values; szfh = sdf['z_form_upper'].values
gelo = np.clip(gzf-gzfl, 0.0, None); gehi = np.clip(gzfh-gzf, 0.0, None)
selo = np.clip(szf-szfl, 0.0, None); sehi = np.clip(szfh-szf, 0.0, None)

OUT = 'figures/'
import os; os.makedirs(OUT, exist_ok=True)

# ============================================================
# FIGURE 1: Formation rate + data overplot (analog of Trenti+15 Fig.3 top)
# ============================================================
fig, ax = plt.subplots(figsize=(10,6))
ax2 = ax.twiny()

ax.fill_between(zp, rn, alpha=0.15, color='gray')
ax.plot(zp, rn, 'k-', lw=2, label='Atomic cooling model (Trenti+15 analog, Planck 2018)')
ax.axvspan(8, 15, alpha=0.08, color='blue', label='Peak formation window (z=8–15)')

y_gems = 0.12; y_spark = 0.06
ax.errorbar(gzf, np.full_like(gzf, y_gems), xerr=[gelo, gehi],
            fmt='o', color='royalblue', ms=7, lw=1.5, capsize=3,
            label='GEMS GCs (z_obs=9.625, n=19)')
ax.errorbar(szf, np.full_like(szf, y_spark), xerr=[selo, sehi],
            fmt='s', color='crimson', ms=8, lw=1.5, capsize=3,
            label='Sparkler GCs (z_obs=1.4, n=5)')

ax.set_xlabel('Formation redshift z_form', fontsize=13)
ax.set_ylabel('Normalized GC formation rate (a.u.)', fontsize=13)
ax.set_xlim(1, 25)
ax.set_ylim(-0.02, 1.15)
ax.legend(fontsize=10, loc='upper right')
ax.set_title('GC Formation Epochs vs. Atomic Cooling Model\n(analog of Trenti+15 Fig. 3 top panel, LCDM Planck 2018)', fontsize=12)

lb_ticks = [13.8, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 6.0, 4.0]
z_for_lb = [float(interpolate.interp1d(lbp, zp, bounds_error=False, fill_value='extrapolate')(lt)) for lt in lb_ticks]
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(z_for_lb)
ax2.set_xticklabels([str(lt) for lt in lb_ticks], fontsize=9)
ax2.set_xlabel('Lookback time (Gyr)', fontsize=11)

plt.tight_layout()
fig.savefig(OUT+'figure1_formation_rate.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure 1 saved.')

# ============================================================
# FIGURE 2: Mass-metallicity relation
# ============================================================
fig, ax = plt.subplots(figsize=(8,6))

ax.errorbar(gdf['log10_mass'], gdf['ZH'],
            xerr=None, yerr=[gdf['ZH_err_lo'], gdf['ZH_err_hi']],
            fmt='o', color='royalblue', ms=7, alpha=0.8, capsize=3, label='GEMS (z=9.625)')
ax.errorbar(sdf['log10_mass'], sdf['ZH'],
            xerr=None, yerr=[sdf['ZH_err_lo'], sdf['ZH_err_hi']],
            fmt='s', color='crimson', ms=9, alpha=0.9, capsize=3, label='Sparkler (z=1.4)')

xg = np.linspace(gdf['log10_mass'].min()-0.1, gdf['log10_mass'].max()+0.1, 100)
slg, icg, rvg, pvg, _ = stats.linregress(gdf['log10_mass'], gdf['ZH'])
ax.plot(xg, slg*xg+icg, '--', color='royalblue', alpha=0.7, lw=1.5,
        label=f'GEMS fit: slope={slg:.2f}, r={rvg:.2f}')
xs = np.linspace(sdf['log10_mass'].min()-0.1, sdf['log10_mass'].max()+0.1, 100)
sls, ics, rvs, pvs, _ = stats.linregress(sdf['log10_mass'], sdf['ZH'])
ax.plot(xs, sls*xs+ics, '--', color='crimson', alpha=0.7, lw=1.5,
        label=f'Sparkler fit: slope={sls:.2f}, r={rvs:.2f}')

allm = np.concatenate([gdf['log10_mass'].values, sdf['log10_mass'].values])
allzh = np.concatenate([gdf['ZH'].values, sdf['ZH'].values])
slc, icc, rvc, pvc, _ = stats.linregress(allm, allzh)
xc = np.linspace(allm.min()-0.1, allm.max()+0.1, 100)
ax.plot(xc, slc*xc+icc, '-', color='gray', alpha=0.5, lw=2,
        label=f'Combined fit: slope={slc:.2f}, r={rvc:.2f}')

ax.set_xlabel('log$_{10}$(M/M$_\odot$)', fontsize=13)
ax.set_ylabel('[Z/H]', fontsize=13)
ax.legend(fontsize=10)
ax.set_title('Mass-Metallicity Relation: GEMS and Sparkler GCs', fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(OUT+'figure2_mass_metallicity.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure 2 saved.')

# ============================================================
# FIGURE 3: CDFs + age distributions + ZH vs z_form
# ============================================================
fig = plt.figure(figsize=(15,5))
gs = gridspec.GridSpec(1, 3, figure=fig)

# Panel a: CDFs
ax1 = fig.add_subplot(gs[0])
sg, cg = np.sort(gzf), np.arange(1,len(gzf)+1)/len(gzf)
ss, cs = np.sort(szf), np.arange(1,len(szf)+1)/len(szf)
ax1.step(sg, cg, where='post', color='royalblue', lw=2, label='GEMS')
ax1.step(ss, cs, where='post', color='crimson', lw=2, label='Sparkler')
ax1.plot(zp, cdfp, 'k--', lw=1.5, alpha=0.7, label='Theory (cumulative)')
ax1.set_xlabel('z_form', fontsize=12); ax1.set_ylabel('CDF', fontsize=12)
ax1.set_title('(a) Cumulative z_form distributions', fontsize=11)
ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3)

# Panel b: Age histograms
ax2 = fig.add_subplot(gs[1])
ax2.hist(gdf['age_Gyr'], bins=8, color='royalblue', alpha=0.6, label='GEMS', density=True)
ax2.hist(sdf['age_Gyr'], bins=5, color='crimson', alpha=0.6, label='Sparkler', density=True)
ax2.set_xlabel('Observed age (Gyr)', fontsize=12); ax2.set_ylabel('Density', fontsize=12)
ax2.set_title('(b) Age distributions', fontsize=11)
ax2.legend(fontsize=9); ax2.grid(True, alpha=0.3)

# Panel c: ZH vs z_form
ax3 = fig.add_subplot(gs[2])
ax3.errorbar(gzf, gdf['ZH'], yerr=[gdf['ZH_err_lo'], gdf['ZH_err_hi']],
             xerr=[gelo, gehi], fmt='o', color='royalblue', ms=6, capsize=2, alpha=0.8, label='GEMS')
ax3.errorbar(szf, sdf['ZH'], yerr=[sdf['ZH_err_lo'], sdf['ZH_err_hi']],
             xerr=[selo, sehi], fmt='s', color='crimson', ms=8, capsize=2, alpha=0.9, label='Sparkler')
allzf2 = np.concatenate([gzf, szf])
allzh2 = np.concatenate([gdf['ZH'].values, sdf['ZH'].values])
slcz, iccz, rvcz, pvcz, _ = stats.linregress(allzf2, allzh2)
xz = np.linspace(allzf2.min()-0.5, allzf2.max()+0.5, 100)
ax3.plot(xz, slcz*xz+iccz, '-', color='gray', alpha=0.6, lw=2,
         label=f'Combined fit: slope={slcz:.3f}')
rho_c, p_c = stats.spearmanr(allzf2, allzh2)
ax3.set_xlabel('z_form', fontsize=12); ax3.set_ylabel('[Z/H]', fontsize=12)
ax3.set_title(f'(c) Metallicity vs formation redshift\nSpearman ρ={rho_c:.2f}, p={p_c:.3f}', fontsize=11)
ax3.legend(fontsize=9); ax3.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(OUT+'figure3_distributions.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure 3 saved.')

# ============================================================
# FIGURE 4: Survival bias — GEMS vs Sparkler mass distributions
# ============================================================
fig, ax = plt.subplots(figsize=(9,6))

ax.hist(sd['gems_log10_masses'], bins=12, color='royalblue', alpha=0.4, density=True, label='GEMS (observed, n=19)')
ax.hist(sd['sparkler_log10_masses'], bins=5, color='crimson', alpha=0.6, density=True, label='Sparkler (observed, n=5)')

alphas = [0.5, 1.0, 2.0]
colors = ['#88c0d0', '#5e81ac', '#2e3440']
for a, c in zip(alphas, colors):
    key = 'alpha_' + str(a).replace('.','p')
    if key in sd:
        ax.hist(sd[key], bins=20, color=c, alpha=0.35, density=True, histtype='step', lw=2,
                label=f'Simulated surviving GEMS (α={a})')

ax.set_xlabel('log$_{10}$(M/M$_\odot$)', fontsize=13)
ax.set_ylabel('Density', fontsize=12)
ax.set_title('Survival Bias: Can Sparkler GCs be the surviving tail of GEMS-like population?', fontsize=12)
ax.legend(fontsize=10); ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(OUT+'figure4_survival_bias.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure 4 saved.')

print('\nAll figures saved to', OUT)
