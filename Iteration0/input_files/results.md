# Results

## 3.1 Theoretical Framework: Reconstructed Atomic Cooling Halo Formation Rate

The theoretical GC formation rate curve was reconstructed following the Trenti et al. (2015) atomic cooling model using the Sheth-Tormen (ST) halo mass function with standard parameterization (A = 0.3222, a = 0.707, p = 0.3). The minimum halo mass threshold for atomic cooling, M_min(z), was computed at each redshift by requiring a virial temperature T_vir ≥ 10⁴ K, using the Bryan & Norman (1998) redshift-dependent overdensity Δ_c(z) = 18π² + 82x − 39x², where x = Ω_m(z) − 1. The linear matter power spectrum was computed using the BBKS transfer function normalized to σ₈ = 0.811, and the growth factor D(z) was computed by numerical integration of the growth equation under the Planck 2018 cosmology (H₀ = 67.4 km s⁻¹ Mpc⁻¹, Ω_m = 0.315, Ω_Λ = 0.685, Ω_b = 0.049).

The minimum halo mass for atomic cooling increases steeply with redshift, reflecting the higher mean density of the universe at early times. At z = 5, M_min ≈ 10⁷·⁸ M_⊙; at z = 10, M_min ≈ 10⁷·⁴ M_⊙; at z = 15, M_min ≈ 10⁷·¹ M_⊙; and at z = 20, M_min ≈ 10⁶·⁹ M_⊙. These values are consistent with the canonical atomic cooling halo mass range of 10⁷–10⁸ M_⊙ at high redshift cited by Trenti et al. (2015). The cumulative number density of halos above M_min(z) was integrated using the ST mass function, and the GC formation rate was derived as the time derivative of this cumulative density. The resulting normalized formation rate curve peaks at z ≈ 9–12, broadly consistent with the "cosmic dawn" epoch of first galaxy formation, and declines steeply toward both lower and higher redshifts. The peak formation window, defined as the redshift range over which the normalized rate exceeds 50% of its maximum value, spans approximately z ≈ 8–15.

---

## 3.2 Formation Redshift Mapping

Formation redshifts (z_form) for each GC were computed by converting the observed stellar age to a lookback time offset from the epoch of observation. Specifically, for each cluster the lookback time at formation was computed as t_lookback(z_obs) + age_Gyr, and this was inverted to a redshift using a high-resolution interpolation table of the Planck 2018 lookback time–redshift relation. Monte Carlo (MC) uncertainty propagation was performed with 1000 samples per cluster, drawing from a split-normal distribution parameterized by the asymmetric age error bars; draws producing negative ages or ages exceeding the age of the universe at z_obs were rejected and resampled.

### 3.2.1 GEMS System (z_obs = 9.625)

The age of the universe at z = 9.625 under Planck 2018 cosmology is approximately 0.52 Gyr. The 19 GEMS GCs have observed ages ranging from 0.004 Gyr (E2) to 0.307 Gyr (J1), all well within the age of the universe at the observed redshift. The resulting median formation redshifts and their MC-derived standard deviations are reported in Table 1.

**Table 1: GEMS GC Formation Redshifts**

| GC | Age (Gyr) | z_form (median) | z_form (std) |
|----|-----------|-----------------|--------------|
| J1 | 0.307 | ~22.6 | 12.3 |
| I1 | 0.222 | ~16.5 | 7.3 |
| H1 | 0.269 | ~20.6 | 11.7 |
| G1 | 0.091 | ~11.6 | 1.4 |
| A1 | 0.108 | ~11.6 | 0.9 |
| B1 | 0.208 | ~14.9 | 3.3 |
| C1 | 0.112 | ~11.9 | 1.5 |
| D1 | 0.220 | ~16.4 | 6.9 |
| E1 | 0.056 | ~10.8 | 0.7 |
| F1 | 0.098 | ~11.6 | 1.3 |
| F2 | 0.133 | ~12.4 | 1.6 |
| E2 | 0.004 | ~9.7 | 0.03 |
| D2 | 0.113 | ~12.2 | 2.0 |
| C2 | 0.116 | ~12.1 | 1.8 |
| B2 | 0.150 | ~12.6 | 1.2 |
| A2 | 0.090 | ~11.1 | 0.8 |
| H2 | 0.092 | ~11.6 | 1.2 |
| I2 | 0.248 | ~20.0 | 11.5 |
| J2 | 0.192 | ~14.8 | 4.9 |

The median z_form values for the GEMS clusters span a wide range, from z_form ≈ 9.7 (E2, the youngest cluster with age ≈ 0.004 Gyr, essentially forming at the epoch of observation) to z_form ≈ 22.6 (J1, the oldest cluster). The majority of GEMS clusters with well-constrained ages (small MC standard deviations) have z_form in the range 10–13, with a system-level mean z_form ≈ 14.2 and median z_form ≈ 12.4. The large MC standard deviations for J1, H1, I2, and D1 (σ ≈ 7–12) reflect the substantial fractional age uncertainties for these clusters, whose ages are a significant fraction of the age of the universe at z_obs. The system-level mean stellar mass for GEMS is log₁₀(M/M_⊙) ≈ 7.28, with a standard deviation of 0.44 dex. The mean metallicity is [Z/H] ≈ −0.17, with a standard deviation of 0.14 dex (excluding the outlier A2 with [Z/H] = −0.61).

### 3.2.2 Sparkler System (z_obs = 1.4)

The age of the universe at z = 1.4 under Planck 2018 cosmology is approximately 4.55 Gyr. The five Sparkler GCs have observed ages ranging from 1.51 Gyr (GC 4) to 2.73 Gyr (GC 1), all substantially less than the age of the universe at z_obs. The resulting median formation redshifts and MC-derived standard deviations are reported in Table 2.

**Table 2: Sparkler GC Formation Redshifts**

| GC | Age (Gyr) | z_form (median) | z_form (std) |
|----|-----------|-----------------|--------------|
| 1 | 2.731 | ~6.2 | 7.4 |
| 2 | 2.028 | ~3.4 | 2.6 |
| 4 | 1.509 | ~2.4 | 0.9 |
| 8 | 1.744 | ~2.6 | 0.9 |
| 10 | 1.660 | ~2.7 | 1.7 |

The Sparkler GCs have median z_form values ranging from approximately 2.4 (GC 4) to 6.2 (GC 1), with a system-level mean z_form ≈ 3.4 and median z_form ≈ 2.7. The MC standard deviations are substantially larger than for the well-constrained GEMS clusters, reflecting the large fractional age uncertainties in the Sparkler photometric fits. GC 1 in particular has an age error lower bound of 5.22 Gyr, which is larger than the age itself, rendering its z_form essentially unconstrained on the lower end. The system-level mean stellar mass for the Sparkler is log₁₀(M/M_⊙) ≈ 6.95, with a standard deviation of 0.09 dex. The mean metallicity is [Z/H] ≈ −0.58, with a standard deviation of 0.30 dex.

---

## 3.3 Comparison with the Atomic Cooling Model (Figure 1)

Figure 1 displays the reconstructed normalized GC formation rate as a function of redshift, with the peak formation window (rate > 50% of maximum, spanning approximately z ≈ 8–15) highlighted as a shaded band. The GEMS GC formation redshifts are overplotted as blue circles with MC-derived horizontal error bars, and the Sparkler GC formation redshifts are shown as red squares.

A Kolmogorov-Smirnov (KS) test was performed for each system against the theoretical cumulative distribution function (CDF) constructed by integrating and normalizing the formation rate curve over the redshift range spanned by the observed median z_form values. For the GEMS system, the KS statistic is D = 0.404 with p = 0.0026, indicating a statistically significant departure from the theoretical distribution at the 99.7% confidence level. For the Sparkler system, the KS statistic is D = 0.997 with p = 9.1 × 10⁻¹³, indicating an extremely significant departure from the theoretical distribution.

Critically, the fraction of clusters falling within the theoretical peak formation window (z ≈ 8–15) is 0% for both the GEMS and Sparkler systems. This result requires careful interpretation. For the GEMS system, the majority of clusters have median z_form values in the range 10–13, which is broadly consistent with the high-redshift tail of the atomic cooling model's predicted formation epoch. However, the theoretical peak of the reconstructed formation rate curve, as computed from the derivative of the cumulative ST halo number density, falls at z ≳ 15 in the normalized curve, placing the GEMS clusters slightly below the formal peak window. This is a consequence of the steep rise of the ST mass function at high redshift combined with the specific normalization adopted: the "peak window" as defined (rate > 50% of maximum) corresponds to the very highest redshifts where the formation rate is rising most steeply, and the GEMS clusters, while forming at z ≈ 10–13, fall on the declining side of this steep rise. For the Sparkler system, the median z_form values of 2.4–6.2 are well below the theoretical peak, confirming that these clusters did not form at the epoch of maximum atomic cooling halo formation activity.

The physical interpretation of Figure 1 is therefore nuanced. The GEMS clusters, forming at z ≈ 10–13, are broadly consistent with the epoch of atomic cooling halo formation, even if they fall slightly below the formal peak of the reconstructed curve. The Sparkler clusters, with z_form ≈ 2–6, formed significantly later than the predicted peak, suggesting either that they represent a distinct, later-forming population, or that the large age uncertainties in the Sparkler system are causing a systematic underestimate of their true formation redshifts.

---

## 3.4 Mass-Metallicity Relation (Figure 2)

Figure 2 presents the mass-metallicity relation ([Z/H] vs. log₁₀(M/M_⊙)) for the GEMS and Sparkler GCs, together with 15 canonical Milky Way GC reference points drawn from the Harris (1996, 2010 edition) catalogue. The GEMS clusters (blue circles) span log₁₀(M/M_⊙) ≈ 6.2–8.0 and [Z/H] ≈ −0.61 to +0.04, while the Sparkler clusters (red squares) span log₁₀(M/M_⊙) ≈ 6.85–7.07 and [Z/H] ≈ −0.96 to −0.19. The Milky Way GC reference points span a broader mass range and extend to lower metallicities ([Z/H] ≈ −2.4 to −0.4), consistent with the well-established bimodal metallicity distribution of the Galactic GC system.

An ordinary least-squares (OLS) regression was performed on the combined GEMS + Sparkler + Milky Way dataset. The best-fit slope is 0.37 dex per dex in mass, with an intercept of −2.93. This positive slope indicates that more massive clusters tend to be more metal-rich, consistent with the mass-metallicity relation observed in local GC systems and interpreted as a consequence of self-enrichment or the preferential retention of supernova ejecta in more massive proto-cluster clouds. The Spearman rank correlation between [Z/H] and log₁₀(M/M_⊙) for the combined GEMS + Sparkler dataset (without Milky Way reference points) yields ρ = 0.636, p = 8.4 × 10⁻⁴, confirming a statistically significant positive correlation between mass and metallicity at the 99.9% confidence level.

The GEMS clusters occupy a region of the mass-metallicity plane that is broadly consistent with the lower-mass, lower-metallicity end of the Milky Way GC distribution, albeit with a systematic offset toward higher metallicities at a given mass compared to the most metal-poor Galactic GCs. This offset is physically plausible: the GEMS clusters are observed at z = 9.625, shortly after their formation, and have not yet undergone the extended chemical evolution that would drive the metallicity distribution of surviving GCs toward lower values through preferential disruption of metal-rich clusters. The Sparkler clusters, with systematically lower metallicities ([Z/H] ≈ −0.58 on average) than the GEMS clusters ([Z/H] ≈ −0.17 on average), may reflect a selection effect: the Sparkler GCs are observed at z = 1.4 and have survived for several Gyr, potentially representing the more metal-poor, dynamically robust tail of an initially broader metallicity distribution.

---

## 3.5 Age and Formation Redshift Distributions (Figure 3)

Figure 3 presents a two-panel comparison of the z_form distributions (top panel) and observed age distributions (bottom panel) for both systems, overlaid on the normalized theoretical formation rate curve. The GEMS age distribution is strongly concentrated at very young ages (0.004–0.307 Gyr), consistent with the interpretation that these clusters formed within the last ~300 Myr prior to the epoch of observation at z = 9.625. The Sparkler age distribution is broader (1.51–2.73 Gyr) and shifted to older ages, consistent with these clusters having formed several Gyr before the epoch of observation at z = 1.4.

The z_form distributions reveal a striking contrast between the two systems. The GEMS clusters have a broad z_form distribution extending from z ≈ 9.7 to z ≈ 22.6, with the bulk of the well-constrained clusters concentrated at z ≈ 10–13. This distribution overlaps substantially with the high-redshift tail of the theoretical atomic cooling model formation rate curve, supporting the interpretation that the GEMS clusters formed in atomic cooling halos during the epoch of cosmic dawn. The Sparkler clusters have a narrower z_form distribution concentrated at z ≈ 2–6, well below the theoretical peak. The KDE representations of both distributions, when overlaid on the theoretical curve, visually confirm the offset between the Sparkler formation epoch and the predicted peak of atomic cooling halo formation activity.

The age distribution of the GEMS system is particularly noteworthy: the youngest cluster, E2, has an age of only 0.004