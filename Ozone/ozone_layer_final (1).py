#!/usr/bin/env python3
"""
ozone_layer.py  --  Sunlight through an OZONE layer lying above an N2 layer
===========================================================================
Run:   python ozone_layer.py                 (needs numpy, matplotlib; data in ./data)
       python ozone_layer.py --du 300 --n2 1e18 --mu0 1.0

THE MODEL
---------
    Sun (blackbody, 5772 K, diluted by distance)
         |
         v        ozone layer  (top)   : column N_O3  [molecules/cm^2]
         v        N2 layer     (below) : column N_N2
         v
    transmitted sunlight

PHYSICS (every step is one line of code below)
----------------------------------------------
1. Source.  The Sun is treated as a blackbody. Planck's law gives the specific intensity
   B_lambda(T); the Sun subtends a small solid angle Omega = pi (R_sun/d)^2, so the irradiance
   at Earth is   F0(lambda) = B_lambda(T_sun) * Omega        [W m^-2 nm^-1]   (integrates to ~1361 W/m^2).
2. Absorption.  A molecule has an effective "target area" for a photon, the cross-section
   sigma(lambda) [cm^2].  Through a column of N molecules/cm^2 the beam is attenuated by
   Beer-Lambert:      F = F_in * exp(-tau / mu0),     tau = N * sigma(lambda).
   mu0 = cos(solar zenith angle): a tilted beam crosses 1/mu0 times more gas.
3. Two layers in series.  Optical depths add, so transmissions multiply:
        F_after_ozone = F0 * exp(-tau_O3/mu0);    F_out = F_after_ozone * exp(-tau_N2/mu0).
   (Hence the ORDER of the layers does not change F_out; it changes the spectrum *between* them.)
4. No scattering, no emission.  Rayleigh scattering is left out on purpose (kept simple); thermal
   emission of cold gas is utterly negligible at these wavelengths.   -> pure absorption.

WHAT THE OZONE ABSORPTION SPECTRUM (the "profile") MEANS
--------------------------------------------------------
A photon is absorbed only if its energy E = hc/lambda matches a jump the molecule can make; a peak
means that jump is likely.  Position = energy match, height = transition probability.
 * HARTLEY band  (~200-310 nm, peak ~255 nm = 4.86 eV, sigma ~1.2e-17 cm^2): by far the strongest.
   A broad continuum (with weak diffuse ripples near the top), nearly temperature independent
   (~2 % over 193-293 K, Serdyuchenko et al. 2014).  Textbook picture (NOT checked in the papers
   used here): the upper state is repulsive, the molecule dissociates (mainly O(1D)+O2), and the
   peak sits where the photon energy matches the upper-state energy at the molecule's equilibrium
   geometry; its width reflects the spread of that geometry (reflection principle).
 * HUGGINS band  (~310-350 nm): ~1e-3..1e-4 of the Hartley strength, with vibrational ripples,
   and STRONGLY temperature dependent: Serdyuchenko et al. attribute temperature dependence to the
   redistribution of ground-state rotational-vibrational populations (hot bands) -- warmer gas has
   more molecules already vibrating, which changes absorption at the band edge.
 * CHAPPUIS band (~400-700 nm, peak ~602 nm = 2.06 eV, sigma ~5e-21 cm^2): ~2000x weaker than
   the Hartley peak, broad, ~1 % temperature dependence; absorbs orange-red light (ozone is
   faintly blue).  The extremely low strength means a low transition probability.
 * N2 (Huffman 1963, 79.8-100 nm): sharp peaks and gaps = electronic-vibrational band structure of
   bound upper states (the b, c4', c5 ... Rydberg/valence systems), unlike ozone's smooth hump.

DATA
----
 N2 : Huffman et al. (1963), 295 K, 79.8-100 nm (139 sparse, irregular points) -> interpolated in
      log(sigma); this is only a coarse ENVELOPE of the real band structure.
 O3 : Serdyuchenko et al. (2014), MPI-Mainz Spectral Atlas, 213.33-1100 nm, 0.01 nm steps, at
      193, 213 and 263 K (tiny negative noise values at long wavelengths are clipped to zero).

ASSUMPTIONS / LIMITS (please keep in the report)
-----------------------------------------------
 A1 Each layer is isothermal at the temperature of its dataset; the three ozone temperatures are
    run as separate scenarios (nothing is interpolated in T).  N2 is fixed at 295 K.
 A2 One uniform layer per gas, fixed column density (no vertical profile).
 A3 Sun = blackbody.  Reasonable in the visible/near-IR, poor in the UV, and wrong in the extreme
    UV (real solar EUV is emission lines) -> the N2 result at 79.8-100 nm is only illustrative.
 A4 Absorption is computed ONLY where data exist: N2 79.8-100 nm, O3 213.33-1100 nm.  Elsewhere in
    50-1150 nm (50-79.8, 100-213.33, 1100-1150 nm) sigma is UNDEFINED: those wavelengths are shaded
    grey, are not modelled (in reality N2 absorbs below 79.6 nm and O2/O3 absorb at 100-213 nm) and
    are excluded from every absorbed-fraction number.
 A5 No scattering (Rayleigh) and no emission.
"""
import os, re, glob, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ user parameters
N_O3_DU = 300.0          # ozone column in Dobson units (1 DU = 2.687e16 molecules/cm^2)
N_N2 = 1.0e18            # N2 column [molecules/cm^2]
MU0 = 1.0                # cos(solar zenith angle); 1.0 = Sun overhead
LAM_MIN, LAM_MAX, DLAM = 50.0, 1150.0, 0.01     # wavelength range and step [nm]
DU = 2.687e16

# ------------------------------------------------------------------ constants (SI)
h, c, kB, sigma_SB = 6.62607015e-34, 2.99792458e8, 1.380649e-23, 5.670374419e-8
R_SUN, AU, T_SUN = 6.957e8, 1.495978707e11, 5772.0
HC_EV_NM = 1239.84198    # hc in eV nm  -> E[eV] = 1239.84 / lambda[nm]


# ------------------------------------------------------------------ source
def planck_nm(lam_nm, T):
    """Planck specific intensity per nm  [W m^-2 sr^-1 nm^-1]."""
    lam = np.asarray(lam_nm, float) * 1e-9
    with np.errstate(over="ignore", divide="ignore"):
        B = 2 * h * c ** 2 / lam ** 5 / np.expm1(h * c / (lam * kB * T))
    return B * 1e-9


def solar_irradiance(lam_nm, T=T_SUN):
    """Irradiance at 1 AU, perpendicular to the beam [W m^-2 nm^-1] = B * Omega, Omega = pi (R/d)^2."""
    return planck_nm(lam_nm, T) * np.pi * (R_SUN / AU) ** 2


# ------------------------------------------------------------------ data
def load_two_columns(path):
    rows = []
    for line in open(path, errors="ignore"):
        p = line.split()
        if len(p) >= 2:
            try:
                rows.append((float(p[0]), float(p[1])))
            except ValueError:
                pass            # header / comment line
    a = np.array(rows)
    return a[np.argsort(a[:, 0])]


def load_data(data_dir):
    n2f = glob.glob(os.path.join(data_dir, "*N2*Huffman*"))
    o3f = glob.glob(os.path.join(data_dir, "*O3*Serdyuchenko*"))
    if not n2f or not o3f:
        raise SystemExit(f"Data files not found in '{data_dir}' (need N2*Huffman* and O3*Serdyuchenko*K*).")
    n2 = load_two_columns(n2f[0])
    o3 = {}
    for f in o3f:
        T = int(re.search(r"(\d{3})K", os.path.basename(f)).group(1))
        a = load_two_columns(f)
        a[:, 1] = np.clip(a[:, 1], 0.0, None)      # remove tiny negative noise
        o3[T] = a
    return n2, dict(sorted(o3.items()))


def sigma_n2(lam, n2):
    """N2 cross-section [cm^2]; log-interpolated between the sparse points; NaN outside the data."""
    s = 10 ** np.interp(lam, n2[:, 0], np.log10(n2[:, 1]))
    return np.where((lam >= n2[0, 0]) & (lam <= n2[-1, 0]), s, np.nan)


def sigma_o3(lam, a):
    """O3 cross-section [cm^2]; linear interpolation of the dense data; NaN outside the data."""
    return np.interp(lam, a[:, 0], a[:, 1], left=np.nan, right=np.nan)


def tau(N, sigma):
    """Optical depth N*sigma; wavelengths without data give 0 (and are masked out of all statistics)."""
    return np.where(np.isfinite(sigma), N * sigma, 0.0)


def integrate(lam, f, lo, hi):
    m = (lam >= lo) & (lam <= hi)
    return np.trapezoid(f[m], lam[m])


# ------------------------------------------------------------------ main
def main(a):
    N_O3 = a.du * DU
    lam = np.arange(LAM_MIN, LAM_MAX + DLAM / 2, DLAM)
    n2, o3 = load_data(a.data)
    temps = list(o3)
    F0 = solar_irradiance(lam)
    s_n2 = sigma_n2(lam, n2)
    s_o3 = {T: sigma_o3(lam, o3[T]) for T in temps}
    has_n2 = np.isfinite(s_n2)
    has_o3 = np.isfinite(s_o3[temps[0]])
    covered = has_n2 | has_o3
    t_n2 = tau(a.n2, s_n2)
    t_o3 = {T: tau(N_O3, s_o3[T]) for T in temps}
    F_o3 = {T: F0 * np.exp(-t_o3[T] / a.mu0) for T in temps}                 # after ozone layer
    F_out = {T: F_o3[T] * np.exp(-t_n2 / a.mu0) for T in temps}              # after N2 layer

    print(f"Ozone column {a.du:g} DU = {N_O3:.3e} cm^-2 | N2 column {a.n2:.2e} cm^-2 | mu0 = {a.mu0:g}")
    # ----- how much of the Sun lies where we have data
    tot = sigma_SB * T_SUN ** 4 * (R_SUN / AU) ** 2
    P_cov = integrate(lam, F0 * covered, LAM_MIN, LAM_MAX)
    P_rng = integrate(lam, F0, LAM_MIN, LAM_MAX)
    print(f"\nSolar constant (blackbody) = {tot:.1f} W/m^2; power in 50-1150 nm = {P_rng:.1f} W/m^2 "
          f"({100*P_rng/tot:.1f} % of total)")
    print(f"Wavelengths WITH data hold {P_cov:.1f} W/m^2 = {100*P_cov/P_rng:.1f} % of the 50-1150 nm power "
          f"(N2 window: {integrate(lam, F0*has_n2, LAM_MIN, LAM_MAX):.2e} W/m^2)")

    # ----- peaks of the ozone spectrum
    print("\nOZONE ABSORPTION BANDS (peak positions measured from the data)")
    print(f"{'T[K]':>5} {'Hartley peak':>14} {'E[eV]':>6} {'sigma':>10} {'FWHM[nm]':>9} {'Chappuis peak':>14} {'sigma':>10}   sigma(330nm)")
    peaks = {}
    for T in temps:
        d = o3[T]; l, s = d[:, 0], d[:, 1]
        mh = (l >= 213.33) & (l <= 320); i = np.argmax(s[mh]); lh, sh = l[mh][i], s[mh][i]
        half = l[mh][s[mh] >= sh / 2]
        ms = np.convolve(s, np.ones(101) / 101, mode="same")          # 1-nm running mean (the band is broad)
        mc = (l >= 450) & (l <= 750); j = np.argmax(ms[mc]); lc, sc = l[mc][j], ms[mc][j]
        s330 = np.interp(330.0, l, s)
        peaks[T] = (lh, sh, lc, sc)
        print(f"{T:5d} {lh:11.2f} nm {HC_EV_NM/lh:6.3f} {sh:10.3e} {half.max()-half.min():9.1f} {lc:11.1f} nm {sc:10.2e}   {s330:.2e}")
    print(f"Hartley/Chappuis peak ratio ({temps[0]} K): {peaks[temps[0]][1]/peaks[temps[0]][3]:.0f}")
    print(f"tau(255 nm) for {a.du:g} DU = {N_O3*peaks[temps[0]][1]:.1f}   -> transmission e^-tau = {np.exp(-N_O3*peaks[temps[0]][1]/a.mu0):.1e}")

    # ----- UV/visible budget for each temperature
    bands = [("UV-C   213.33-280", 213.33, 280), ("UV-B   280-315", 280, 315), ("UV-A   315-400", 315, 400),
             ("visible 400-700", 400, 700), ("near-IR 700-1100", 700, 1100)]
    print("\nPOWER BUDGET inside the ozone data window (W/m^2 perpendicular to the beam)")
    print(f"{'band':<18} {'incident':>9} " + " ".join(f"{'T='+str(T)+'K out':>12} {'absorbed':>9}" for T in temps))
    for name, lo, hi in bands:
        inc = integrate(lam, F0, lo, hi)
        row = f"{name:<18} {inc:9.2f} "
        for T in temps:
            out = integrate(lam, F_out[T], lo, hi)
            row += f"{out:12.3f} {100*(1-out/inc):8.2f}% "
        print(row)
    inc_all = integrate(lam, F0 * has_o3, LAM_MIN, LAM_MAX)
    row = f"{'whole O3 window':<18} {inc_all:9.2f} "
    for T in temps:
        out = integrate(lam, F_out[T] * has_o3, LAM_MIN, LAM_MAX)
        row += f"{out:12.3f} {100*(1-out/inc_all):8.2f}% "
    print(row)
    print(f"   = {100*(inc_all-integrate(lam,F_out[temps[1]]*has_o3,LAM_MIN,LAM_MAX))/tot:.2f} % of the TOTAL solar constant is absorbed "
          f"({temps[1]} K case)")
    inc_n2 = integrate(lam, F0 * has_n2, LAM_MIN, LAM_MAX); out_n2 = integrate(lam, F_out[temps[0]] * has_n2, LAM_MIN, LAM_MAX)
    print(f"N2 window 79.8-100 nm: incident {inc_n2:.3e} W/m^2, transmitted {out_n2:.3e} W/m^2 "
          f"(mean tau = {np.nanmean(t_n2[has_n2]):.1f}, mean transmission {np.mean(np.exp(-t_n2[has_n2]/a.mu0)):.2e})")

    # ----- built-in checks
    print("\nCHECKS")
    ok = []
    def check(name, good, detail=""):
        ok.append(good); print(f"  {'PASS' if good else 'FAIL'}  {name}  {detail}")
    lg = np.logspace(0, 6, 400001); Fsb = np.pi * planck_nm(lg, T_SUN) * 1.0
    check("Stefan-Boltzmann: integral of pi*B dlambda = sigma T^4", abs(np.trapezoid(Fsb, lg) / (sigma_SB * T_SUN ** 4) - 1) < 2e-3,
          f"(ratio {np.trapezoid(Fsb, lg) / (sigma_SB * T_SUN ** 4):.5f})")
    lm = np.linspace(200, 2000, 200001); lmax = lm[np.argmax(planck_nm(lm, T_SUN))]
    check("Wien: lambda_max * T = 2.8978e-3 m K", abs(lmax * 1e-9 * T_SUN / 2.897772e-3 - 1) < 1e-3, f"(lambda_max = {lmax:.1f} nm)")
    chk = F_out[temps[0]] * 1.0; direct = F0 * np.exp(-(t_o3[temps[0]] + t_n2) / a.mu0)
    check("Two layers multiply: F_out = F0 exp(-(tau_O3+tau_N2)/mu0)", np.allclose(chk, direct, rtol=1e-12, atol=0))
    A1 = integrate(lam, F0 * (1 - np.exp(-tau(1e12, s_o3[temps[0]]))), 213.33, 1100)
    A2 = integrate(lam, F0 * (1 - np.exp(-tau(2e12, s_o3[temps[0]]))), 213.33, 1100)
    check("Thin-layer limit: absorbed power doubles when a tiny column doubles", abs(A2 / A1 - 2) < 1e-3, f"(ratio {A2/A1:.5f})")
    check("Hartley peak at 255 +- 0.5 nm", all(abs(peaks[T][0] - 255.4) < 0.5 for T in temps))
    check("Chappuis peak at 602 +- 3 nm", all(abs(peaks[T][2] - 602) < 3 for T in temps))
    check("Warmer gas absorbs more at 330 nm (Huggins hot-band effect)", np.interp(330, o3[temps[-1]][:, 0], o3[temps[-1]][:, 1]) > np.interp(330, o3[temps[0]][:, 0], o3[temps[0]][:, 1]))
    check("Hartley peak changes < 3 % between coldest and warmest file", abs(peaks[temps[-1]][1] / peaks[temps[0]][1] - 1) < 0.03)
    check("No negative cross-sections after clipping", all((o3[T][:, 1] >= 0).all() for T in temps) and (n2[:, 1] > 0).all())
    print(f"\n{sum(ok)}/{len(ok)} checks passed")

    # ------------------------------------------------------------ figures
    os.makedirs(a.out, exist_ok=True)
    BLUE, ORANGE, TEAL, GRAY, RED = "#2a6fdb", "#e8710a", "#12a594", "#9aa0a6", "#d62728"
    cols = dict(zip(temps, (BLUE, TEAL, ORANGE)))
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.18, "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False})

    def shade_nodata(ax):
        for lo, hi in ((LAM_MIN, 79.8), (100, 213.33), (1100, LAM_MAX)):
            ax.axvspan(lo, hi, color=GRAY, alpha=0.18, lw=0)

    # Fig 1: cross-sections with the bands explained
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    shade_nodata(ax)
    ax.plot(n2[:, 0], n2[:, 1], "o", color=RED, ms=3, label="N$_2$ (Huffman 1963, 295 K)")
    for T in temps:
        d = o3[T]; m = d[:, 1] > 0
        ax.plot(d[m, 0], d[m, 1], color=cols[T], lw=0.8, label=f"O$_3$ {T} K (Serdyuchenko 2014)")
    lh, sh, lc, sc = peaks[temps[0]]
    ax.annotate(f"HARTLEY peak\n{lh:.1f} nm = {HC_EV_NM/lh:.2f} eV\nσ = {sh:.2e} cm²\nbroad continuum\n(dissociation)", (lh, sh), (140, 2e-14), fontsize=8.5, arrowprops=dict(arrowstyle="->", color="k"))
    ax.annotate("HUGGINS band\nweak, rippled,\nstrongly T-dependent", (325, 1.5e-20), (350, 3e-18), fontsize=8.5, arrowprops=dict(arrowstyle="->", color="k"))
    ax.annotate(f"CHAPPUIS peak\n{lc:.0f} nm = {HC_EV_NM/lc:.2f} eV\nσ = {sc:.1e} cm²\n(~{sh/sc:.0f}× weaker)", (lc, sc), (720, 3e-18), fontsize=8.5, arrowprops=dict(arrowstyle="->", color="k"))
    ax.text(88, 2e-22, "N$_2$ bands:\nsharp peaks/gaps", color=RED, fontsize=8.5, ha="center")
    ax.text(150, 1e-23, "no data\n(not modelled)", color="gray", fontsize=8.5, ha="center")
    ax.text(62, 1e-23, "no data", color="gray", fontsize=8, ha="center"); ax.text(1125, 1e-23, "no\ndata", color="gray", fontsize=8, ha="center")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(LAM_MIN, LAM_MAX); ax.set_ylim(1e-24, 1e-10)
    ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Cross-section σ (cm$^2$/molecule)")
    ax.set_title("Absorption profiles of N$_2$ and O$_3$: where and why the peaks are"); ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=8.5, ncol=4)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig1_cross_sections.png"), dpi=170); plt.close(fig)

    # Fig 2: temperature effect
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
    for T in temps:
        d = o3[T]; m = (d[:, 0] > 245) & (d[:, 0] < 266)
        ax[0].plot(d[m, 0], d[m, 1] * 1e17, color=cols[T], lw=1, label=f"{T} K")
    ax[0].set_xlabel("Wavelength (nm)"); ax[0].set_ylabel("σ (10$^{-17}$ cm$^2$)"); ax[0].set_title("(a) Hartley peak: almost no T dependence"); ax[0].legend()
    base = temps[0]; w = np.ones(51) / 51
    for T in temps[1:]:
        d = o3[T]; l = d[:, 0]; sm = np.convolve(d[:, 1], w, mode="same"); sb = np.convolve(o3[base][:, 1], w, mode="same")
        m = (l > 295) & (l < 350) & (sb > 0)
        ax[1].plot(l[m], sm[m] / sb[m], color=cols[T], label=f"σ({T} K)/σ({base} K)")
    ax[1].axhline(1, color="k", lw=0.6); ax[1].set_ylim(0.8, 2.2); ax[1].set_xlabel("Wavelength (nm)")
    ax[1].set_ylabel("Ratio to the coldest file"); ax[1].set_title("(b) Huggins band: strong T dependence"); ax[1].legend()
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig2_temperature_effect.png"), dpi=170); plt.close(fig)

    # Fig 3: sunlight spectra
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    m = (lam >= 213.33) & (lam <= 1100)
    ax[0].plot(lam[m], F0[m], color=GRAY, label="top of atmosphere")
    ax[0].plot(lam[m], F_out[temps[1]][m], color=BLUE, label=f"after O$_3$ + N$_2$ ({temps[1]} K)")
    ax[0].set_xlabel("Wavelength (nm)"); ax[0].set_ylabel("Irradiance (W m$^{-2}$ nm$^{-1}$)"); ax[0].legend(); ax[0].set_title("(a) Whole ozone window")
    m = (lam >= 213.33) & (lam <= 400)
    ax[1].semilogy(lam[m], F0[m], color=GRAY, label="top of atmosphere")
    for T in temps:
        ax[1].semilogy(lam[m], np.maximum(F_out[T][m], 1e-30), color=cols[T], lw=1, label=f"{T} K")
    ax[1].set_ylim(1e-10, 2); ax[1].set_xlabel("Wavelength (nm)"); ax[1].set_title(f"(b) UV zoom, {a.du:g} DU: Hartley blocks UV-C/B"); ax[1].legend(fontsize=8)
    m = (lam >= 79.8) & (lam <= 100)
    ax[2].semilogy(lam[m], F0[m], color=GRAY, label="incident"); ax[2].semilogy(lam[m], np.maximum(F_out[temps[0]][m], 1e-30), color=RED, lw=0.8, label="after N$_2$ layer (floored at 1e-30)"); ax[2].set_ylim(1e-30, 10)
    ax[2].set_xlabel("Wavelength (nm)"); ax[2].set_title(f"(c) N$_2$ window (N = {a.n2:.0e} cm$^{{-2}}$)"); ax[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig3_sunlight_spectra.png"), dpi=170); plt.close(fig)
    print(f"\nFigures written to {a.out}/")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--data", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"))
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "output"))
    ap.add_argument("--du", type=float, default=N_O3_DU, help="ozone column [Dobson units]")
    ap.add_argument("--n2", type=float, default=N_N2, help="N2 column [cm^-2]")
    ap.add_argument("--mu0", type=float, default=MU0, help="cos(solar zenith angle)")
    main(ap.parse_args())
