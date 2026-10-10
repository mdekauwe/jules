"""Offline size check: linear (current l_som_rhizo_series) vs matric-flux-
potential rhizosphere at FR-Pue, from a best-fit run's half-hourly output.
Units: head h (m, suction > 0), K (kg m-2 s-1 per unit head gradient),
B (m-1 ground), flux kg m-2 s-1. Mirrors smc_ext soil_to_root_conductance
and JULES vG (hh_from_sthu, hyd_con_vg)."""
import sys
import numpy as np
import netCDF4 as nc
from scipy.optimize import brentq

run = sys.argv[1]
dz = np.array([0.0500, 0.0711, 0.1010, 0.1435, 0.2039, 0.2898])
b, sathh, ks, sm_sat = 4.97, 0.816, 0.016, 0.30
rootd, lma, lai_bal = 0.8592, 0.226, 1.82
r_r, rho_r = 0.29e-3, 0.31e3
RHOG = 1000.0 * 9.8
M_H2O = 0.018

z2 = np.cumsum(dz); z1 = z2 - dz
f_root = (np.exp(-z1 / rootd) - np.exp(-z2 / rootd)) / (1 - np.exp(-z2[-1] / rootd))
m_root = lma * lai_bal
B = 4 * np.pi / np.log(rho_r * dz / (f_root * m_root)) * f_root * m_root / (rho_r * np.pi * r_r**2)
rho = (f_root * m_root / (rho_r * dz)) ** -0.5          # r_s / r_r


def se_of_h(h):
    h = np.maximum(h, 0.0)
    return (1.0 + (h / sathh) ** ((b + 1) / b)) ** (-1.0 / (b + 1))


def k_of_se(s):
    s = np.clip(s, 0.0, 1.0)
    sd = np.clip(s, 0.05, 0.95)
    k = ks * sd**0.5 * (1 - (1 - sd**(b + 1)) ** (1 / (b + 1))) ** 2
    lo = s < 0.05
    k = np.where(lo, k / 0.05 * s, k)
    hi = s > 0.95
    k = np.where(hi, k + (ks - k) / 0.05 * (s - 0.95), k)
    return k


def h_of_se(s):
    s = np.maximum(s, 0.01)
    return sathh * (s ** (-b - 1) - 1.0) ** (b / (b + 1))


# Phi(h) = int_h^hmax K dh (hmax = h at sthu_min), tabulated on log h
hmax = h_of_se(0.01)
hg = np.concatenate([[0.0], np.logspace(-4, np.log10(hmax), 20000)])
kg = k_of_se(se_of_h(hg))
cum = np.concatenate([[0.0], np.cumsum(0.5 * (kg[1:] + kg[:-1]) * np.diff(hg))])
phig = cum[-1] - cum


def phi(h):
    return np.interp(h, hg, phig)


def solve(hs, E, hr_on):
    """Root suction h_r with sum B (Phi(hs) - Phi(h_r)) = E."""
    def f(hr):
        d = phi(hs) - phi(hr)
        if not hr_on:
            d = np.maximum(d, 0.0)
        return np.sum(B * d) - E
    lo = 0.0 if hr_on else hs.min()
    if f(hmax) < 0:
        return np.nan                               # beyond soil supply
    return brentq(f, lo, hmax, xtol=1e-6)


ds = nc.Dataset(f"{run}/outputs/run.out.nc")
t = nc.num2date(ds["time"][:], ds["time"].units, only_use_cftime_datetimes=False)
sm = ds["SoilMoist"][:, :, 0, 0]
tv = ds["TVeg"][:, 0, 0]
prz = ds["psi_root_zone_pft"][:, 0, 0, 0]
pl = ds["psi_leaf_pft"][:, 0, 0, 0]
kp = ds["kplant_pft"][:, 0, 0, 0]
to_mmol_mpa = 1e6 / RHOG / M_H2O * 1e3            # kg m-2 s-1 m-1 -> mmol m-2 s-1 MPa-1

print("layers: f_root", np.round(f_root, 3), " rho=r_s/r_r", np.round(rho, 1),
      " ln rho", np.round(np.log(rho), 2))
print("steady-rate/steady-state geometry factor ratio:",
      np.round(np.log(rho) / (rho**2 * np.log(rho) / (rho**2 - 1) - 0.5), 3))

hdr = ("period", "n", "psi_s*", "psi_r lin", "psi_r MFP", "psi_r MFP noHR",
       "k_lin", "k'_MFP", "k'noHR", "kplant", "psi_leaf", "Esupply/E", "psi_rz(model)", "E mm/h")
print("\nmedians over 11-14 UTC half-hours with TVeg > 0.01 mm/h; psi MPa, k mmol m-2 s-1 MPa-1")
print(" | ".join(hdr))
periods = [("2003 JJA", 2003, (6, 7, 8)), ("2006 JJA", 2006, (6, 7, 8)),
           ("2007 JJA", 2007, (6, 7, 8)), ("2005 Aug", 2005, (8,)), ("2012 Aug", 2012, (8,)),
           ("2004 Apr-May (wet)", 2004, (4, 5)), ("2010 Apr-May (wet)", 2010, (4, 5))]
for name, yr, months in periods:
    rows = []
    for i in range(len(t)):
        ti = t[i]
        if ti.year != yr or ti.month not in months or not 11 <= ti.hour < 14:
            continue
        E = float(tv[i])
        if E < 0.01 / 3600:
            continue
        s = np.clip(sm[i] / (1000 * dz * sm_sat), 0.01, 1.0)
        hs = h_of_se(s)
        k = B * k_of_se(s)
        ksum = k.sum()
        h_rz = np.sum(k * hs) / ksum
        hr_lin = h_rz + E / ksum
        hr_mfp = solve(hs, E, True)
        hr_nohr = solve(hs, E, False)
        kp_mfp = np.sum(B * k_of_se(se_of_h(hr_mfp))) if np.isfinite(hr_mfp) else 0.0
        kp_nohr = (np.sum(B * k_of_se(se_of_h(hr_nohr)) * (hs <= hr_nohr))
                   if np.isfinite(hr_nohr) else 0.0)
        esup = np.sum(B * phi(hs)) / E
        rows.append([-h_rz * RHOG / 1e6, -hr_lin * RHOG / 1e6, -hr_mfp * RHOG / 1e6,
                     -hr_nohr * RHOG / 1e6, ksum * to_mmol_mpa, kp_mfp * to_mmol_mpa,
                     kp_nohr * to_mmol_mpa, float(kp[i]), float(pl[i]), esup, float(prz[i]), E*3600])
    r = np.array(rows)
    if len(r) == 0:
        print(name, 0); continue
    med = np.nanmedian(r, axis=0)
    print(f"{name} | {len(r)} | " + " | ".join(f"{v:.3g}" for v in med))
    # driest 5 % of steps by psi_s*
    q = r[r[:, 0] <= np.nanpercentile(r[:, 0], 5)]
    print(f"   driest 5% | {len(q)} | " + " | ".join(f"{v:.3g}" for v in np.nanmedian(q, axis=0)))
