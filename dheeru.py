"""
1D Radiative Transfer Model: Solar Blackbody -> Constant-Density N2 Slab
=========================================================================

Physical picture
-----------------
A blackbody at the Sun's photospheric temperature (T = 5778 K) sits right
next to a 1D slab of pure molecular nitrogen (N2) at constant number
density N, with path length L. We propagate the blackbody spectrum
through the slab using the Beer-Lambert law:

        I(lambda) = I0(lambda) * exp( -sigma(lambda) * N * L )

I0(lambda)      : incident blackbody spectral radiance (Planck's law)
sigma(lambda)   : N2 photoabsorption cross-section (wavelength-dependent)
N               : N2 number density (constant, cm^-3)
L               : slab path length (cm)

Why N2 shows NO discrete lines in the IR/visible
-------------------------------------------------
N2 is homonuclear -> no permanent dipole moment and no dipole change on
vibration -> all rotational and vibrational electric-dipole transitions
are forbidden. So in the IR/visible/near-UV, N2 is essentially
transparent (no discrete structure at all).

Where N2 DOES show real discrete structure
--------------------------------------------
In the vacuum UV (VUV), N2 has strong, fully dipole-ALLOWED electronic
transitions:
  * Worley-Jenkins / Rydberg-valence bands (b, b', c'_4 <- X), ~80-100 nm
    -> the strongest, sharpest "picket fence" of absorption lines,
       several overlapping Rydberg series converging toward the
       ionization limit.
  * Lyman-Birge-Hopfield (LBH) bands (a <- X), ~127-170 nm
    -> weaker (electric-dipole forbidden, g<-g), sparser vibronic bands.
  * Photoionization continuum below ~79.6 nm (IP = 15.58 eV)
    -> smooth, non-discrete absorption edge (bound-free, not bound-bound).

A 5778 K blackbody has just enough flux in this VUV window for these
transitions to matter, which is why this project uses the solar
temperature rather than Earth's 288 K surface temperature.

IMPORTANT CAVEAT ABOUT THE CROSS-SECTION MODEL USED HERE
-----------------------------------------------------------
Real, lab-measured N2 VUV cross-sections (the kind used in actual
aeronomy/thermosphere models) come from tabulated databases such as
Fennelly & Torr (1992) or the MPI-Mainz UV/VIS Spectral Atlas. This
script does NOT have network access to fetch those tables, so instead
it builds a PHYSICALLY-MOTIVATED SYNTHETIC MODEL of the N2 VUV
cross-section:
  - Three overlapping Rydberg series (standing in for the b, b', c'_4
    states) built from the Rydberg formula E_n = IP - Ry/(n-delta)^2,
    correctly converging toward the real N2 ionization limit.
  - A separate, weaker, roughly evenly-spaced vibrational progression
    standing in for the LBH (a<-X) system.
  - A smooth power-law photoionization continuum below threshold.

This reproduces the CORRECT QUALITATIVE PHYSICS (line positions
converging correctly, relative strengths, discrete-vs-continuum
transition at the real ionization energy) and is suitable for an
educational/illustrative simulation. For a research-grade or
publication-grade result, replace `n2_cross_section()` with an
interpolation over a real tabulated dataset.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Output folder: figures are saved next to this script, in a
# subfolder called "output" -- this works the same way on Windows,
# macOS, and Linux (no hardcoded /home/... or C:\... paths).
# ----------------------------------------------------------------------
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # __file__ isn't defined in some interactive environments (e.g. a
    # Jupyter cell); fall back to the current working directory.
    SCRIPT_DIR = os.getcwd()

OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ----------------------------------------------------------------------
# Physical constants (SI)
# ----------------------------------------------------------------------
h  = 6.62607015e-34      # Planck constant, J s
c  = 2.99792458e8        # speed of light, m/s
kB = 1.380649e-23        # Boltzmann constant, J/K

T_SUN = 5778.0           # solar photosphere temperature, K

EV_TO_NM = 1239.841984   # hc in eV*nm  ->  wavelength[nm] = EV_TO_NM / E[eV]


def planck_lambda(wavelength_m, T):
    """Spectral radiance B_lambda(T) [W / (m^2 sr m)] via Planck's law."""
    x = h * c / (wavelength_m * kB * T)
    # guard against overflow for very short wavelengths
    x = np.clip(x, None, 700)
    return (2.0 * h * c**2) / (wavelength_m**5) / (np.expm1(x))


def nm_to_eV(wl_nm):
    return EV_TO_NM / wl_nm


def eV_to_nm(E_eV):
    return EV_TO_NM / E_eV


# ----------------------------------------------------------------------
# Synthetic N2 VUV photoabsorption cross-section model
# ----------------------------------------------------------------------
IP = 15.581     # N2 ionization potential, eV  (-> 79.58 nm threshold)
RY = 13.6057    # Rydberg constant, eV


def rydberg_line_centers(n_values, quantum_defect, IP=IP):
    """Line energies (eV) for a Rydberg series converging to IP."""
    n_values = np.asarray(n_values, dtype=float)
    return IP - RY / (n_values - quantum_defect) ** 2


def lorentzian(x, x0, width):
    """Normalized-ish Lorentzian profile (peak = 1)."""
    return (width / 2) ** 2 / ((x - x0) ** 2 + (width / 2) ** 2)


def build_line_list():
    """
    Returns a list of (center_nm, width_nm, peak_sigma_cm2) tuples
    describing every discrete N2 VUV line in the model.
    """
    lines = []

    # --- Rydberg-valence (Worley-Jenkins) series: b, b', c'_4 <- X ---
    # Three overlapping series with different quantum defects, standing
    # in for the three mixed excited states. n runs from 3 up to ~9;
    # higher members crowd together near the ionization limit, exactly
    # as in the real spectrum.
    series_params = [
        dict(defect=1.15, n_max=9, peak_sigma=3.0e-17, base_width_nm=0.06,  label="b<-X series"),
        dict(defect=0.95, n_max=9, peak_sigma=4.5e-17, base_width_nm=0.05,  label="b'<-X series"),
        dict(defect=0.75, n_max=9, peak_sigma=6.0e-17, base_width_nm=0.045, label="c'4<-X series"),
    ]
    for sp in series_params:
        n_vals = np.arange(3, sp["n_max"] + 1)
        energies_eV = rydberg_line_centers(n_vals, sp["defect"])
        centers_nm = eV_to_nm(energies_eV)
        for i, wl0 in enumerate(centers_nm):
            if wl0 < IP_THRESH_NM:   # only keep lines below the ionization edge (bound states)
                continue
            # oscillator strength / cross-section falls off for higher n,
            # and lines get slightly broader near the series limit
            # (increasing predissociation / level density)
            decay = np.exp(-0.35 * i)
            sigma_peak = sp["peak_sigma"] * decay
            width = sp["base_width_nm"] * (1 + 0.4 * i)
            lines.append((wl0, width, sigma_peak))

    # --- Lyman-Birge-Hopfield (a <- X) vibrational progression ---
    # Symmetry-forbidden -> much weaker; roughly evenly spaced in energy
    # (vibrational progression rather than a converging Rydberg series).
    lbh_v = np.arange(0, 9)
    E0, dE = 7.35, 0.255   # eV: band origin and vibrational spacing (approx.)
    for v in lbh_v:
        E = E0 + dE * v
        wl0 = eV_to_nm(E)
        sigma_peak = 4.0e-19 * np.exp(-0.12 * v)   # weak, slowly varying envelope
        width = 0.35
        lines.append((wl0, width, sigma_peak))

    return lines


IP_THRESH_NM = eV_to_nm(IP)   # ~79.58 nm ionization threshold

LINE_LIST = build_line_list()


def n2_cross_section(wl_nm):
    """
    Total N2 photoabsorption cross-section [cm^2] at wavelength(s) wl_nm.
    = sum of discrete Lorentzian lines (bound-bound)
      + smooth power-law continuum below the ionization threshold (bound-free)
    """
    sigma = np.zeros_like(wl_nm, dtype=float)

    # discrete lines
    for wl0, width, peak in LINE_LIST:
        sigma += peak * lorentzian(wl_nm, wl0, width)

    # smooth photoionization continuum for wl < threshold (shorter wavelength)
    continuum_peak = 2.2e-17  # cm^2, typical order of magnitude at threshold
    below = wl_nm < IP_THRESH_NM
    # simple smooth rise-then-decay continuum shape referenced to threshold
    sigma[below] += continuum_peak * (wl_nm[below] / IP_THRESH_NM) ** 3

    return sigma


# ----------------------------------------------------------------------
# Build the incident blackbody spectrum
# ----------------------------------------------------------------------
# Broad context view (visible -> VUV) to show where N2 is transparent
wl_broad_nm = np.linspace(50, 1000, 4000)
I0_broad = planck_lambda(wl_broad_nm * 1e-9, T_SUN)

# Fine VUV grid (this is where all the discrete N2 structure lives)
wl_vuv_nm = np.linspace(50, 200, 300_000)   # fine spacing to resolve narrow lines
I0_vuv = planck_lambda(wl_vuv_nm * 1e-9, T_SUN)

sigma_vuv = n2_cross_section(wl_vuv_nm)     # cm^2

# ----------------------------------------------------------------------
# Slab parameters -> Beer-Lambert transmission
# ----------------------------------------------------------------------
N_density = 5.0e12     # N2 number density, cm^-3  (representative thermosphere-ish value)
L_path    = 2.0e6       # path length, cm  (20 km, to make line saturation visible)

tau_vuv = sigma_vuv * N_density * L_path
transmittance = np.exp(-tau_vuv)
I_transmitted_vuv = I0_vuv * transmittance

# Apply the same (essentially zero) absorption to the broad-context curve
sigma_broad = n2_cross_section(wl_broad_nm)
tau_broad = sigma_broad * N_density * L_path
I_transmitted_broad = I0_broad * np.exp(-tau_broad)

# ----------------------------------------------------------------------
# PLOT 1: Broad context - incident vs transmitted, visible -> VUV
# ----------------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(10, 5.5))
ax1.plot(wl_broad_nm, I0_broad, color="darkorange", lw=1.8, label="Incident (5778 K blackbody)")
ax1.plot(wl_broad_nm, I_transmitted_broad, color="firebrick", lw=1.2, label="Transmitted through N$_2$ slab")
ax1.set_yscale("log")
ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel(r"Spectral radiance $B_\lambda$ (W m$^{-2}$ sr$^{-1}$ m$^{-1}$)")
ax1.set_title("Solar blackbody vs. transmitted spectrum through N$_2$ slab\n(visible down to VUV - note: flat/unattenuated except deep in the UV)")
ax1.axvspan(80, 100, color="royalblue", alpha=0.12, label="Rydberg-valence (b,b',c') region")
ax1.axvspan(127, 170, color="seagreen", alpha=0.12, label="LBH band region")
ax1.legend(loc="lower right", fontsize=9)
ax1.grid(alpha=0.3, which="both")
fig1.tight_layout()
fig1.savefig(os.path.join(OUTPUT_DIR, "plot1_broad_context.png"), dpi=150)

# ----------------------------------------------------------------------
# PLOT 2: N2 cross-section vs wavelength (the "raw" discrete spectrum)
# ----------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(10, 5.5))
ax2.plot(wl_vuv_nm, sigma_vuv, color="navy", lw=0.8)
ax2.axvline(IP_THRESH_NM, color="black", ls="--", lw=1,
            label=f"Ionization threshold ({IP_THRESH_NM:.1f} nm)")
ax2.set_xlabel("Wavelength (nm)")
ax2.set_ylabel(r"N$_2$ photoabsorption cross-section $\sigma$ (cm$^2$)")
ax2.set_title("Synthetic N$_2$ VUV absorption cross-section\n(discrete Rydberg-valence + LBH lines, converging toward ionization limit)")
ax2.legend(loc="upper right", fontsize=9)
ax2.grid(alpha=0.3)
fig2.tight_layout()
fig2.savefig(os.path.join(OUTPUT_DIR, "plot2_cross_section.png"), dpi=150)

# ----------------------------------------------------------------------
# PLOT 3: Zoomed transmitted spectrum - the actual "dips" result
# ----------------------------------------------------------------------
fig3, (ax3a, ax3b) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

ax3a.plot(wl_vuv_nm, I0_vuv, color="darkorange", lw=1, label="Incident I$_0(\\lambda)$")
ax3a.plot(wl_vuv_nm, I_transmitted_vuv, color="firebrick", lw=1, label="Transmitted I($\\lambda$)")
ax3a.set_ylabel(r"Spectral radiance")
ax3a.set_title("N$_2$ slab transmission in the VUV: discrete absorption lines (dips)")
ax3a.legend(loc="upper left", fontsize=9)
ax3a.grid(alpha=0.3)

ax3b.plot(wl_vuv_nm, transmittance, color="teal", lw=0.8)
ax3b.axvline(IP_THRESH_NM, color="black", ls="--", lw=1, label="Ionization threshold")
ax3b.set_xlabel("Wavelength (nm)")
ax3b.set_ylabel("Transmittance  I/I$_0$")
ax3b.set_ylim(0, 1.05)
ax3b.legend(loc="lower left", fontsize=9)
ax3b.grid(alpha=0.3)

fig3.tight_layout()
fig3.savefig(os.path.join(OUTPUT_DIR, "plot3_transmission_dips.png"), dpi=150)

print("Done. Figures saved in:", OUTPUT_DIR)
print(" -", os.path.join(OUTPUT_DIR, "plot1_broad_context.png"))
print(" -", os.path.join(OUTPUT_DIR, "plot2_cross_section.png"))
print(" -", os.path.join(OUTPUT_DIR, "plot3_transmission_dips.png"))
print(f"Ionization threshold wavelength: {IP_THRESH_NM:.3f} nm")
print(f"Number of discrete lines in model: {len(LINE_LIST)}")