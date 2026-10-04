# Xylem impairment (memory model) and recovery: review, 2026-10-04

Branch `xylem_impairment` at f720dab. Runs are in
`runs/roses/FR_Pue/FR_Pue_dev_slow_recovery/best_seg` (off, rec3, rec5;
root_psi_crit = psi_close = -6 MPa, root-zone driver, segments, flat Ci +
Newton). Offline emulation: `best_seg/recovery_forms_offline.py` and
`plots/recovery_forms_offline.png`.

V = verified (code read or numbers computed); I = inferred.

## 1. Code correctness

No bug changes the current results. These are the confirmed issues:

1. **Linear recovery restarts its clock on any new damage.** This is a
   design flaw. (V) `xylem_impairment_memory_mod.f90:767`:
   `ximpair_plc_dam = 1 - kcap/kmax` whenever kcap falls. PLC_dam is the
   total current loss, not the increment, so the rate becomes
   (current loss)/R. Example: at 3 % residual loss, a new damage to 3.1 %
   cuts the rate from 19/3 to 3.1/3 %/yr and adds about 3 years of tail.
   Every event therefore gets a full R-year tail whatever its size, and
   recovery hits 0 with a corner (2013 and 2022 in the plot).
2. **Stale comment.** (V) `xylem_impairment_memory_mod.f90:598-601` says
   "LAI ... not held in the dump". `ximpair_lai_prev` and all the other
   state *are* dumped and read (`write_dump_mod.F90:390-417`,
   `populate_var.inc:2034`, `required_vars_for_configuration_mod.F90:290-294`).
3. **The growth clock's magnitude is unphysical.** It only works because
   it is normalised. (V) The NPP term renews 2.4 live-sapwood pools per
   year at FR-Pue (0.3 NPP / (eta_sl canht LAI) = 0.3 NPP / about
   0.11 kg C m-2). The LAI term adds 0.43 per year. In the slow scheme
   only the relative timing matters, so this is harmless. In the fast
   scheme (rec_years = 0) the magnitude *is* the recovery rate, which is
   why that scheme was too fast.
4. **The damage memory is not spun up.** (V) The update returns in
   spin-up (`xylem_impairment_mod.f90`, `update_xylem_impairment`,
   IF(is_spinup)). The main run starts intact, and the running mean starts
   at the main run with a bias correction (correct, V). This is minor here:
   the first event is mid-2000 at about 1 %.

Checked and found OK (V):
- **Sentinel.** kmax_impaired = 0 means intact, and is mapped to kmax_pft
  at `init_vars_tmp_mod.F90:336`.
- **Basis.** kcap is on the leaf basis (kmax_pft). In the two-leaf path
  kmax_sun/shd = k_max * lai * share and kmax_ref = kmax_pft * lai * share,
  so `kcap_frac = kmax/kmax_ref = kcap/kmax_pft` (`stom_opt_jls_mod.f90:1231`).
- **Segments.** The cap applies to segments 2-3 only
  (`xylem_hydraulics_cumulative_weibull_jls_mod.f90:925`), with a
  two-piece integral that is continuous at psi_cap.
- **Order of operations per step:** recovery, refill, reset, damage, then
  the kcrit floor (`:761`). One update per step, after the two-leaf
  iteration (`sf_stom_jls_mod.F90:3279`).
- **Damage points.** Damage acts only where apar > 0 (daytime). With the
  root driver this makes no difference.
- **Timestep.** The running mean uses a = dt/(R yr), so it is independent
  of the timestep.
- **Driver.** `psi_root` is the raw root-zone psi, not the gravity-adjusted
  psi_src (`:616`). That is a 0.05 MPa difference, negligible.

### Why impairment slightly increases GPP (V)

It is not the change of cost reference. Before the first damage
(2000-01-01 to 2000-07-14) impairment is already on and the cost
reference is already the intact zero-flow k (`stom_opt_jls_mod.f90:419`).
Yet |dGPP| < 2e-6 g C m-2 d-1 over that period, so the two references
coincide. The effect is soil water carry-over. Damaged plants transpire
slightly less in spring (2018 Apr-Jun: dT -0.0007 to -0.002 mm d-1, dSM
up to +0.25 mm). They spend that water in July-September (2018 Aug:
dGPP +0.011 g C m-2 d-1). The relative effect is 4e-4 at PLC 10-20 %.

On the formulation: the zero-flow intact series k at psi_soil is
Sperry et al. (2017)'s own definition of k_max for the cost, so it is
right. Max-over-samples is only an approximation of it. (I) With a
damage cap one could argue the cost should be zero above psi_cap, since
no new embolism occurs there. The current scheme, which uses the intact
k at the impaired psi, still charges a cost there. That is a defensible
"total loss" view and makes a negligible difference at FR-Pue.

## 2. Exponential vs linear: quantitative comparison

All forms run on the same growth clock. g is the growth in the step in
units of a typical year, i.e. f_new / (running mean x 1 yr), with
bias correction:

- **lin (current):** PLC <- max(PLC - PLC_dam g/R, 0). Fully recovered
  after R typical years. The mean residence of the loss is R/2.
- **exp, growth clock:** PLC <- PLC exp(-g/tau).
  - tau = R/2 gives the same integrated legacy per unit damage as lin.
  - tau = R/3 gives 95 % recovered after R.
- **exp, time only:** PLC <- PLC exp(-dt/tau), tau = R/2, with no growth.

**Emulation.** Damage was replayed daily from the off run's root-zone psi
on the stem Weibull curve (P50 -6.66, P88 -7.75 MPa), as a cap that only
falls, floored at kcrit. Validation of the linear form against JULES:
rec3 RMSE 0.23 pts, bias -0.17, r 0.999; rec5 RMSE 0.25 pts, r 0.999 (V).

| form | R | mean PLC 2018-21 (%) | 2017 peak | t to 50 % of peak (yr) | t to 10 % (yr) | mean 2019-20 | mean 2022-24 |
|---|---|---|---|---|---|---|---|
| lin | 3 | 8.4 | 18.8 | 1.54 | 4.27* | 7.8 | 6.3 |
| exp growth, tau = R/2 | 3 | 7.1 | 18.8 | 1.08 | 4.47* | 6.5 | 5.7 |
| exp growth, tau = R/3 | 3 | 5.8 | 18.8 | 0.65 | 3.87* | 5.2 | 4.9 |
| exp time, tau = R/2 | 3 | 7.0 | 18.8 | 1.04 | 4.37* | 6.4 | 5.8 |
| lin | 5 | 10.7 | 18.8 | 2.54 | 4.53* | 10.8 | 7.2 |
| exp growth, tau = R/2 | 5 | 8.8 | 18.8 | 1.73 | not reached* | 8.2 | 6.7 |
| exp growth, tau = R/3 | 5 | 7.5 | 18.8 | 1.28 | 4.57* | 6.9 | 5.9 |
| exp time, tau = R/2 | 5 | 8.8 | 18.8 | 1.73 | not reached* | 8.2 | 6.8 |
| no recovery | - | 18.8 | 18.8 | - | - | 18.8 | 18.8 |

\* The 10 % times are set by the 2019 (about 9 %) and 2022 (about 11 %)
re-damage, not by the form itself.

Findings (V):
- **Peaks are identical.** Damage dominates them.
- **Forms differ by 1-3 PLC points in the years after an event.** lin
  carries the most legacy for a given R, because it front-loads the loss.
  exp(tau = R/2) and lin give the same integral per event *only without
  re-damage*. With re-damage, lin's clock restarts (section 1.1), so lin
  is higher: 2019-20 is 7.8 vs 6.5 (R = 3) and 10.8 vs 8.2 (R = 5).
- **Back-to-back events.** For 2016 -> 2017, the forms differ by at most
  1 pt (2016-17 mean 6.6-7.7 at R = 3), because 2017 overwrites 2016. For
  2017 -> 2019 and 2022 -> 2023, lin is 1.3-2.6 pts above exp(R/2),
  because of the restart.
- **The growth clock does not matter at FR-Pue.** exp on the growth clock
  and exp on time (same tau) differ by at most 0.1 pt. The interannual CV
  of annual renewal is 0.08, and 85 % of it is the NPP term. Seasonally,
  renewal is 47 % in March-June and 9 % in August-September.
- **R = 3 vs 5 shifts the post-2017 mean by about 2 pts** in all forms.
- **The forms are not separable in the fluxes.** The dGPP from the whole
  impairment is under 0.05 %.

## 3. Paschalis et al. (2023, GCB 30: e17022), the T&C-d scheme

From the PDF, pp. 4-6, 10-12:

- **Damage.** Two nodes, xylem (stem) and leaf, each with a Weibull-type
  vulnerability curve: k_x = exp(-q_x |psi_x|^p_x) k_x^max (Eq. 5) and
  k_l (Eq. 6). The damage drivers are the **xylem-node psi_x for k_x and
  psi_l for k_l**, not psi_soil (p. 5, Sec. 2.6).
- **Memory: an age-cohort structure** (Eq. 14):
  k_x = integral over age a of k_x^a(psi_x^min(a)) p_a da. Each cohort
  remembers the minimum psi it has experienced. There is **no refilling**
  (p. 5). The same scheme applies to k_l.
- **Recovery: only by building new tissue.** The cohort age distribution
  follows a McKendrick-type PDE (Eq. 15) with boundary p_0(t) = lambda_x
  (Eq. 16):
  - New tissue is built at rate lambda_x and is born full.
  - Turnover mu_x(a, t) removes the *oldest* tissue first.
  - Both lambda_x and mu_x come from T&C's dynamic-vegetation carbon pools
    (living sapwood, leaves) (p. 6).
  - So the recovery time is the sapwood (or leaf) renewal and turnover
    time. It is emergent; there is no rec_years parameter, and no NSC
    term.
- **Feedback.** The damaged k_x and k_l enter the soil-xylem-leaf flux
  equations (Eqs. 3-4) and hence the stomatal reduction f_l (Eqs. 8-9).
- **Parameters.** These are in Table S1 (supplement, not seen). The
  vulnerability curves were digitised from published curves.
- **Result at FR-Pue.** T&C-HC-d strongly *overestimated* the 2003 -> 2004
  legacy (Fig. 5e, p. 12; GPP in 2004 about half of observed). Its skill
  fell (Table 2, p. 14: FR-Pue GPP r2 0.73 -> 0.53, KGE 0.74 -> 0.55).
  - The authors attribute this to having no change in carbon allocation
    conditional on xylem damage: the plant does not prioritise restoring
    xylem (p. 12).
  - Damage variants were worse at all sites (p. 15).
  - They cite Page et al. (2023): drought flux legacies are rare beyond
    about 6 months.

**Comparison with ours.**
- Cohort memory with FIFO turnover is a one-pool limit we don't have.
- Our fast branch, `kcap += (kmax - kcap) f_renew` (`:707`), is the
  *well-mixed* version of their renewal: dilution by new tissue. It is an
  exponential on the growth clock with tau = 1/(annual renewal fraction).
- Their scheme ties recovery to the actual renewal rate. This produced
  too much legacy at this very site, because sapwood turnover is slow.
- Ours sets the timescale explicitly (R) and uses growth only as a pacing
  clock. That avoids their failure mode and keeps R calibratable.

**Adopt?**
- No to the cohort PDE: it costs state and parameters and is unobservable
  at FR-Pue.
- Yes to two parts:
  1. The stem-psi driver (theirs is psi_x) together with a separate leaf
     damage driven by psi_l.
  2. The well-mixed renewal form, i.e. exponential on the growth clock,
     which is their scheme's natural one-pool reduction.
- Their FR-Pue result argues *for* short R (2-3 yr or less), consistent
  with Heinrich et al. (2026: no flux legacy after 2017).

## 4. Literature on the recovery form

(I; from memory, not re-fetched here.)
- Anderegg et al. (2015, Science): ring-width legacies of 1-4 years,
  largest in year 1 and decaying. That is a front-loaded shape, closer to
  exponential than linear.
- Kannenberg et al. (2019, 2020): legacies appear in ring width but are
  much weaker or absent in GPP and fluxes. Growth legacies depend on
  carbon allocation, not hydraulics alone.
- Trugman et al. (2018): incomplete hydraulic recovery, with recovery via
  new xylem, can produce multi-year legacies.
- Mackay et al. (2015): TREES caps the curve at the minimum psi seen
  (our cap form).
- Sperry/Venturas (2017): no refilling under tension.

Overall this supports a cap with no refilling, recovery by growth with a
front-loaded (exponential-like) decay, and a timescale of 1-4 years.
**Observability at FR-Pue:** fluxes cannot constrain the form or R. Only
the PLC trajectory or ring width could; Puéchabon has dendrometer and
ring-width data (I).

## 5. Recommendations (ranked)

1. **Exponential recovery on the growth clock, tau = R/2.** Effort: 0.5 day.
   - It removes PLC_dam, and with it the restart artefact (1.1) and the
     corner at zero.
   - It is the one-pool limit of Paschalis's renewal and matches the
     front-loaded literature shape.
   - With tau = R/2 it keeps the integrated legacy of the current lin for
     an isolated event, so the meaning of R carries over.
2. **Damage driver option 4 = stem psi from the segment solver.** Effort:
   1 day.
   - Expose the stem-segment outlet psi (`psi_in` after iseg = 2 in
     `leaf_psi_segments_jls`) of the chosen sample. Then pass it like
     psi_leaf to `update_xylem_impairment`.
   - The root-zone driver ignores the midday drop through the root
     segment, which carries half the resistance, so it under-damages the
     stem.
   - (I) Expect larger PLC in 2017 and more flux effect. This is a
     sensitivity test, not a default yet: the model's root zone is still
     0.8-1.3 MPa too dry.
3. **Fix the stale comment at `:598-601`.** Effort: trivial.
4. **(Optional) A separate leaf damage, reset at leaf flush (LAI term)
   and stem damage recovering with NPP.** Effort: 2 days.
   - This is Paschalis's two-tissue split.
   - It is not identifiable at FR-Pue (evergreen *Q. ilex*, gradual leaf
     turnover). Defer.
5. **Keep refilling (ximpair_tau_rec) off, and no NSC-limited recovery.**
   - There are no data to constrain either.
   - Recovery during low-growth (drought) years is already slowed by the
     NPP term.

## 6. Proposed formulation (option 1)

New namelist `ximpair_rec_form` (jules_vegetation, integer: 1 = linear
(default, current), 2 = exponential). Replace the linear block at
`xylem_impairment_memory_mod.f90:728-732` with:

```
plc(:) = 1.0 - kcap(:) / kmax_pts(:)
SELECT CASE (ximpair_rec_form)
CASE (ximpair_rec_linear)
  plc(:) = MAX(plc(:) - ximpair_plc_dam(:,pft) * growth(:)                    &
                        / ximpair_rec_years(pft), 0.0)
CASE (ximpair_rec_exp)
  ! e-folding of the loss in typical years of growth; tau = R/2 keeps the
  ! integrated legacy of an isolated event equal to the linear form's.
  plc(:) = plc(:) * EXP( -2.0 * growth(:) / ximpair_rec_years(pft) )
END SELECT
kcap(:) = kmax_pts(:) * (1.0 - plc(:))
```

In equations: dPLC/dg = -PLC/tau with tau = R/2. The growth clock stays
as it is (`:713-726`). ximpair_plc_dam stays as a diagnostic only.

Plumbing (same pattern as `l_ximpair_rec_growth`):
`jules_vegetation_mod.F90` declaration near `:545`, the NAMELIST list,
n_int, the my_nml copy in and out, a print, and a check that the value is
in {1, 2}.

## 7. Minimal tests (no more than 4 runs)

1. rec_form = 2, R = 3, root driver. Compare with the emulated "exp growth
   R/2" line (expect RMSE < 0.3 pts) and with lin rec3.
2. rec_form = 2, R = 3, stem-psi driver (option 2). This gives the size
   of the driver effect on PLC and fluxes.
3. off with the same code: confirm it is bit-identical to the current off.
4. (Optional) rec_form = 2, R = 1.5 against Heinrich's no-legacy
   constraint, with the stem driver.

## 8. Note for later: carbon-limited recovery (needs TRIFFID)

Decision (2026-10-04): leave refilling off; keep the fixed-timescale
recovery (`ximpair_rec_years`, linear or exponential via
`ximpair_rec_form`) for now. The natural next step is a carbon-limited
recovery, and that needs TRIFFID (`l_triffid`):

- Today the growth clock (new LAI + `ximpair_wood_alloc` * NPP, or
  TRIFFID wood with `ximpair_growth_basis = 3`) is normalised by its own
  running mean. It only sets the *timing* of recovery; the *amount* per
  year is fixed by `ximpair_rec_years`.
- Carbon-limited version: take the renewed fraction from the new sapwood
  TRIFFID actually allocates (basis 3), without the running-mean
  normalisation, so a year with little carbon for wood recovers little.
  Optionally give the repair an explicit carbon cost taken out of the wood
  allocation (a trade-off between repair and growth).
- Motivation: Paschalis et al. (2023) attribute their overestimated FR-Pue
  2003 -> 2004 legacy to having no extra carbon allocated to xylem repair.
  JULES has no NSC pool, so the TRIFFID wood allocation is the place to
  hook it.
- Prerequisites: a TRIFFID-on FR-Pue setup with sensible wood growth, and
  the growth-clock magnitude issue in section 1 (item 3) fixed: the raw
  NPP term renews about 2.4 sapwood pools a year.
