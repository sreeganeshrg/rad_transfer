"""
==============================================================================
VOIGT LINE PROFILE - a beginner-friendly simulation
==============================================================================

WHAT THIS CODE DOES
--------------------
This program builds up a Voigt spectral line profile step by step:

    Gaussian (Doppler) profile  +  Lorentzian (pressure) profile
                    |
                    v
              Voigt profile
                    |
                    v
   absorption coefficient -> optical depth -> Beer-Lambert transmission

It is written with plain functions, no classes, and no unnecessary
abstraction, so that every line can be read and modified directly.

CONVENTIONS USED IN THIS CODE (read this before changing anything!)
---------------------------------------------------------------------
1.  Everything is in terms of ORDINARY frequency nu (Hz), never angular
    frequency omega = 2*pi*nu (rad/s). If you ever compare this code to a
    textbook that uses omega, you must insert factors of 2*pi.

2.  delta_nu_D ("Doppler width") is defined as the 1/e HALF-width,

        delta_nu_D = (nu0/c) * sqrt(2*kB*T/m)

    This is the convention used by Rybicki & Lightman, "Radiative
    Processes in Astrophysics" (1979/1986), Chapter 10, and is the same
    quantity HITRAN calls the Doppler width. It is NOT the FWHM.
    FWHM_Gaussian = 2*sqrt(ln 2) * delta_nu_D  (about 1.665 * delta_nu_D).

3.  gamma ("Lorentzian width") is defined as the HALF-width at half
    maximum (HWHM) of the Lorentzian, so that the Lorentzian FWHM is
    2*gamma. This matches the HITRAN convention for the pressure
    broadening half-width (see Gordon et al. 2022, HITRAN2020, JQSRT
    277, 1-82).

4.  All three profiles (Gaussian, Lorentzian, Voigt) are normalized so
    that integral over all frequency = 1. They have units of 1/Hz.

5.  The dimensionless Voigt parameters are
        x = (nu - nu0) / delta_nu_D          (dimensionless detuning)
        a = gamma / delta_nu_D               (dimensionless Voigt parameter)

REFERENCES (see the end of the accompanying explanation for full details)
---------------------------------------------------------------------
 - Rybicki, G. B., & Lightman, A. P. (1986). Radiative Processes in
   Astrophysics. Wiley. (Doppler broadening, absorption coefficient,
   optical depth, Beer-Lambert law)
 - Gordon, I. E., et al. (2022). The HITRAN2020 molecular spectroscopic
   database. JQSRT, 277, 1-82. (Voigt line-shape convention, line
   strength S, HWHM convention for gamma)
 - Schreier, F. (1992). The Voigt and complex error function: A
   comparison of computational methods. JQSRT, 48(5), 743-762.
   (relation between the Voigt function and the Faddeeva function)
 - Olivero, J. J., & Longbothum, R. L. (1977). Empirical fits to the
   Voigt line width: A brief review. JQSRT, 17(2), 233-236.
   (approximate Voigt FWHM formula used in Section 16)
 - SciPy documentation for scipy.special.wofz (Faddeeva function)
==============================================================================
"""

# ------------------------------------------------------------------------
# 1. IMPORT LIBRARIES
# ------------------------------------------------------------------------
import numpy as np                     # arrays and numerical operations
import matplotlib.pyplot as plt        # plotting
from scipy.special import wofz         # Faddeeva function, used for the Voigt profile


# ------------------------------------------------------------------------
# 2. PHYSICAL CONSTANTS (SI units)
# ------------------------------------------------------------------------
c_light = 2.99792458e8      # speed of light, m/s (exact, SI definition)
k_B = 1.380649e-23          # Boltzmann constant, J/K (exact, SI definition)
amu = 1.66053906660e-27     # atomic mass unit, kg (CODATA)


# ------------------------------------------------------------------------
# 3. USER-ADJUSTABLE PARAMETERS
# ------------------------------------------------------------------------
# --- physical parameters -------------------------------------------------
nu0 = 2.0e13                 # line-center frequency, Hz
                              # (illustrative value: ~15 micron infrared
                              # transition, similar in scale to the CO2
                              # bending-mode band that matters for Earth's
                              # atmosphere. Not a real, catalogued line.)

temperature = 300.0          # gas temperature, K
mass_amu = 44.0               # molecular mass, atomic mass units
                              # (44 amu = CO2-like mass, illustrative)
mass = mass_amu * amu        # molecular mass, kg

gamma = 2.5e7                 # Lorentzian HWHM, Hz
                              # (manually set; in a real atmosphere this
                              # would come from pressure/collisional
                              # broadening, but here it is just a free
                              # parameter you can change by hand)

# --- radiative-transfer demonstration parameters (Section 15) -----------
line_strength = 1.0e-7        # arbitrary demo units (illustrative only,
                              # NOT a real HITRAN line strength - see
                              # Section 15 for the unit discussion)
column_density = 1.0e15       # arbitrary demo units (illustrative only)

# --- numerical parameters -------------------------------------------------
half_width_in_doppler_widths = 40   # how far the frequency grid extends,
                                     # in units of delta_nu_D, on each
                                     # side of nu0
number_of_points = 4000             # number of points in the frequency grid

# --- plotting parameters ---------------------------------------------------
freq_unit_scale = 1.0e6      # divide frequencies by this to plot in MHz
freq_unit_name = "MHz"       # label for the plot x-axes


# ------------------------------------------------------------------------
# 4. CALCULATE THE DOPPLER WIDTH
# ------------------------------------------------------------------------
# Physics: thermal motion of the absorbing molecules Doppler-shifts the
# frequency each molecule absorbs at. Averaging over the Maxwell-Boltzmann
# velocity distribution turns a single frequency nu0 into a Gaussian
# spread of frequencies with 1/e half-width delta_nu_D.
#
# Reference: Rybicki & Lightman (1986), Ch. 10.
delta_nu_D = (nu0 / c_light) * np.sqrt(2.0 * k_B * temperature / mass)

# Python meaning: this is one line of ordinary NumPy arithmetic. nu0,
# c_light, k_B, temperature and mass are all plain numbers (floats), so
# np.sqrt just computes the square root of a number.
#
# Physics meaning: bigger T -> molecules move faster on average -> bigger
# spread of Doppler shifts -> bigger delta_nu_D. Bigger mass -> molecules
# move MORE SLOWLY at a given temperature (equipartition: (1/2)m<v^2> is
# fixed by T, not m) -> smaller delta_nu_D. That is why light molecules
# (like H2) have broader Doppler profiles than heavy ones (like CO2) at
# the same temperature.
#
# Units check: [nu0/c] = Hz / (m/s) = s/m. [sqrt(k_B*T/m)] = sqrt(J/kg) =
# sqrt(m^2/s^2) = m/s. Multiplying: (s/m)*(m/s) = dimensionless... but we
# want Hz. Let's be careful: nu0 has units Hz = 1/s, c has units m/s, so
# nu0/c has units (1/s)/(m/s) = 1/m. Multiplying by a velocity (m/s) gives
# 1/s = Hz. Correct.


# ------------------------------------------------------------------------
# 5. CALCULATE THE VOIGT PARAMETER a
# ------------------------------------------------------------------------
# a tells you, at a glance, whether Doppler or Lorentzian broadening
# dominates: a << 1 means the line looks almost Gaussian, a >> 1 means it
# looks almost Lorentzian.
a_parameter = gamma / delta_nu_D


# ------------------------------------------------------------------------
# 6. DEFINE THE FREQUENCY GRID
# ------------------------------------------------------------------------
# We build the grid centered on nu0, wide enough (in units of delta_nu_D)
# to show the Lorentzian wings, which fall off much more slowly than the
# Gaussian ones.
half_width_hz = half_width_in_doppler_widths * delta_nu_D
nu_grid = np.linspace(nu0 - half_width_hz, nu0 + half_width_hz, number_of_points)

# Python meaning: np.linspace(start, stop, N) makes an array of N evenly
# spaced numbers from start to stop (inclusive). nu_grid is therefore a
# 1-D NumPy array of frequencies, in Hz.
#
# Physics meaning: this is simply the list of frequencies at which we
# will evaluate the line profile - it plays the same role as the x-axis
# of a spectrum.


# ------------------------------------------------------------------------
# 7. GAUSSIAN (DOPPLER) PROFILE
# ------------------------------------------------------------------------
def gaussian_profile(nu, nu0, delta_nu_D):
    """
    Normalized Gaussian (Doppler) line profile.
    Reference: Rybicki & Lightman (1986), Eq. 10.xx (Doppler profile).

    Inputs:
        nu        - frequency or array of frequencies, Hz
        nu0       - line-center frequency, Hz
        delta_nu_D - Doppler width (1/e half-width), Hz
    Returns:
        phi_G(nu), units of 1/Hz, normalized so that
        integral of phi_G over all nu equals 1.
    """
    prefactor = 1.0 / (delta_nu_D * np.sqrt(np.pi))
    return prefactor * np.exp(-((nu - nu0) / delta_nu_D) ** 2)


# ------------------------------------------------------------------------
# 8. LORENTZIAN PROFILE
# ------------------------------------------------------------------------
def lorentzian_profile(nu, nu0, gamma):
    """
    Normalized Lorentzian line profile (HWHM = gamma).

    Inputs:
        nu    - frequency or array of frequencies, Hz
        nu0   - line-center frequency, Hz
        gamma - Lorentzian half-width at half maximum (HWHM), Hz
    Returns:
        phi_L(nu), units of 1/Hz, normalized so that
        integral of phi_L over all nu equals 1.
    """
    return (gamma / np.pi) / ((nu - nu0) ** 2 + gamma ** 2)


# ------------------------------------------------------------------------
# 9. VOIGT PROFILE (via the Faddeeva function)
# ------------------------------------------------------------------------
def voigt_profile(nu, nu0, delta_nu_D, gamma):
    """
    Normalized Voigt line profile: the convolution of a Gaussian and a
    Lorentzian, evaluated using the Faddeeva function w(z) = exp(-z^2)*erfc(-i*z).

    scipy.special.wofz(z) computes exactly this function w(z) for complex z.
    If we choose z = x + i*a, where x = (nu-nu0)/delta_nu_D and a =
    gamma/delta_nu_D, then it is a standard mathematical result (see
    Schreier 1992, JQSRT 48, 743) that

        Re[w(x + i*a)] = K(x, a)

    where K(x,a) is the (unnormalized) Voigt function

        K(x,a) = (a/pi) * integral_{-inf}^{inf} exp(-y^2) / ((x-y)^2+a^2) dy .

    K(x,a) already has the property that, for fixed a,
        integral over x of K(x,a) dx = sqrt(pi),
    so dividing by (delta_nu_D * sqrt(pi)) gives a profile normalized to 1
    in frequency space. This is why the prefactor below looks just like
    the Gaussian prefactor.

    Why use wofz instead of doing the convolution integral by hand?
    The direct convolution integral of a Gaussian and a Lorentzian has no
    closed form in terms of elementary functions - it can only be written
    as this complex error function. wofz is a fast, accurate, well-tested
    numerical implementation of that special function, so calling it is
    both simpler and far more numerically reliable than integrating the
    convolution directly.

    Inputs:
        nu         - frequency or array of frequencies, Hz
        nu0        - line-center frequency, Hz
        delta_nu_D - Doppler width (1/e half-width), Hz
        gamma      - Lorentzian HWHM, Hz
    Returns:
        phi_V(nu), units of 1/Hz, normalized so that
        integral of phi_V over all nu equals 1 (to numerical precision).
    """
    x = (nu - nu0) / delta_nu_D
    a = gamma / delta_nu_D
    z = x + 1j * a
    return np.real(wofz(z)) / (delta_nu_D * np.sqrt(np.pi))


# Evaluate all three profiles on the frequency grid, using the default
# parameters set in Section 3.
phi_gaussian = gaussian_profile(nu_grid, nu0, delta_nu_D)
phi_lorentzian = lorentzian_profile(nu_grid, nu0, gamma)
phi_voigt = voigt_profile(nu_grid, nu0, delta_nu_D, gamma)


# ------------------------------------------------------------------------
# 10. NORMALIZATION CHECK
# ------------------------------------------------------------------------
# Physically, phi(nu) redistributes the strength of a transition across
# frequency - broadening spreads the same total absorption/emission
# probability over a range of frequencies, it does not create or destroy
# it. That physical statement is exactly the mathematical requirement
# that integral of phi(nu) dnu = 1. If a profile is not normalized to 1,
# any line strength S multiplying it would give the wrong total
# absorption.
#
# IMPORTANT SUBTLETY: the Lorentzian (and therefore the Voigt) wings die
# off only as 1/(nu-nu0)^2, much more slowly than the Gaussian's
# exp(-x^2). This means the plotting grid from Section 6 (only
# +/-40 delta_nu_D wide, chosen so the plots are readable) is NOT wide
# enough to capture all the area in the wings, and integrating phi_voigt
# over that grid alone would under-count the normalization by a percent
# or so. To check normalization properly we build a second, much wider
# grid just for this integral - a good example of how a range that looks
# "wide enough" for a plot can still be too narrow for an accurate
# numerical integral.
integration_half_width_in_doppler_widths = 2000
integration_number_of_points = 40000
nu_grid_for_integration = np.linspace(
    nu0 - integration_half_width_in_doppler_widths * delta_nu_D,
    nu0 + integration_half_width_in_doppler_widths * delta_nu_D,
    integration_number_of_points,
)
phi_voigt_for_integration = voigt_profile(nu_grid_for_integration, nu0, delta_nu_D, gamma)

if hasattr(np, "trapezoid"):
    voigt_integral = np.trapezoid(phi_voigt_for_integration, nu_grid_for_integration)  # NumPy >= 2.0
else:
    voigt_integral = np.trapz(phi_voigt_for_integration, nu_grid_for_integration)      # older NumPy

print(f"Integral of Voigt profile = {voigt_integral:.6f}  (should be close to 1)")


# ------------------------------------------------------------------------
# 11. PLOT 1 - BASIC VOIGT PROFILE
# ------------------------------------------------------------------------
detuning_grid = (nu_grid - nu0) / freq_unit_scale  # frequency detuning, in MHz

plt.figure(figsize=(7, 5))
plt.plot(detuning_grid, phi_voigt, color="tab:blue")
plt.axvline(0, color="gray", linestyle="--", linewidth=1, label="line center")
plt.xlabel(f"Frequency detuning  (nu - nu0)  [{freq_unit_name}]")
plt.ylabel("Normalized Voigt profile  [1/Hz]")
plt.title("Plot 1: Basic Voigt line profile")
plt.legend()
plt.tight_layout()
plt.savefig("plot1_basic_voigt.png", dpi=130)
plt.close()


# ------------------------------------------------------------------------
# 12. PLOT 2 - GAUSSIAN vs LORENTZIAN vs VOIGT
# ------------------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(detuning_grid, phi_gaussian, label="Gaussian (Doppler only)", linestyle="--")
plt.plot(detuning_grid, phi_lorentzian, label="Lorentzian (pressure only)", linestyle=":")
plt.plot(detuning_grid, phi_voigt, label="Voigt (both combined)", linewidth=2)
plt.xlabel(f"Frequency detuning  (nu - nu0)  [{freq_unit_name}]")
plt.ylabel("Normalized profile  [1/Hz]")
plt.title("Plot 2: Gaussian vs Lorentzian vs Voigt")
plt.legend()
plt.tight_layout()
plt.savefig("plot2_gaussian_lorentzian_voigt.png", dpi=130)
plt.close()


# ------------------------------------------------------------------------
# 13. PLOT 3 - EFFECT OF TEMPERATURE
# ------------------------------------------------------------------------
temperature_values = [100.0, 300.0, 600.0]  # K

plt.figure(figsize=(7, 5))
for T_i in temperature_values:
    delta_nu_D_i = (nu0 / c_light) * np.sqrt(2.0 * k_B * T_i / mass)
    phi_i = voigt_profile(nu_grid, nu0, delta_nu_D_i, gamma)
    plt.plot(detuning_grid, phi_i, label=f"T = {T_i:.0f} K")
plt.xlabel(f"Frequency detuning  (nu - nu0)  [{freq_unit_name}]")
plt.ylabel("Normalized Voigt profile  [1/Hz]")
plt.title("Plot 3: Effect of temperature (gamma fixed)")
plt.legend()
plt.tight_layout()
plt.savefig("plot3_temperature_effect.png", dpi=130)
plt.close()


# ------------------------------------------------------------------------
# 14. PLOT 4 - EFFECT OF THE LORENTZIAN WIDTH gamma
# ------------------------------------------------------------------------
gamma_values = [0.2 * delta_nu_D, 1.0 * delta_nu_D, 5.0 * delta_nu_D]  # Hz

plt.figure(figsize=(7, 5))
for gamma_i in gamma_values:
    phi_i = voigt_profile(nu_grid, nu0, delta_nu_D, gamma_i)
    plt.plot(detuning_grid, phi_i, label=f"gamma = {gamma_i/freq_unit_scale:.2f} {freq_unit_name}")
plt.xlabel(f"Frequency detuning  (nu - nu0)  [{freq_unit_name}]")
plt.ylabel("Normalized Voigt profile  [1/Hz]")
plt.title("Plot 4: Effect of the Lorentzian width gamma (T fixed)")
plt.legend()
plt.tight_layout()
plt.savefig("plot4_gamma_effect.png", dpi=130)
plt.close()


# ------------------------------------------------------------------------
# 15. PLOT 5 - EFFECT OF THE DIMENSIONLESS VOIGT PARAMETER a
# ------------------------------------------------------------------------
# Here we directly control a = gamma/delta_nu_D by choosing gamma relative
# to the fixed default delta_nu_D, to show the Gaussian-to-Lorentzian
# transition explicitly.
a_values_demo = [0.05, 0.5, 2.0, 10.0]

plt.figure(figsize=(7, 5))
for a_i in a_values_demo:
    gamma_i = a_i * delta_nu_D
    phi_i = voigt_profile(nu_grid, nu0, delta_nu_D, gamma_i)
    plt.plot(detuning_grid, phi_i, label=f"a = {a_i:g}")
plt.xlabel(f"Frequency detuning  (nu - nu0)  [{freq_unit_name}]")
plt.ylabel("Normalized Voigt profile  [1/Hz]")
plt.title("Plot 5: Effect of the Voigt parameter a = gamma / delta_nu_D")
plt.legend()
plt.tight_layout()
plt.savefig("plot5_voigt_parameter_a.png", dpi=130)
plt.close()


# ------------------------------------------------------------------------
# 16. SIMPLE RADIATION-TRANSFER DEMONSTRATION
# ------------------------------------------------------------------------
# This section is intentionally kept separate and simple - it is only a
# first taste of how the Voigt profile feeds into radiative transfer, not
# a full atmospheric model.
#
# Concepts (see Rybicki & Lightman 1986, Ch. 1, for the general
# absorption-coefficient / optical-depth / Beer-Lambert framework):
#
#   alpha_nu = S * phi_V(nu)          absorption coefficient per unit
#                                      column density, units: length^2
#                                      (a cross section, spread over
#                                      frequency by phi_V)
#
#   tau_nu   = N_column * alpha_nu    optical depth for a slab of column
#                                      density N_column (dimensionless)
#
#   T_nu     = exp(-tau_nu)           Beer-Lambert transmission
#
# UNIT NOTE (please read): a real HITRAN line strength S has units of
# cm/molecule (because HITRAN works in wavenumber, not frequency). Here,
# because our profile phi_V(nu) is normalized in frequency [1/Hz], the
# "line_strength" parameter below must be understood as an integrated
# cross section per molecule with units [length^2 * Hz]. We are treating
# it purely as an illustrative, hand-set number (not a real molecular
# line strength) so you can see how S, column density, and phi_V combine
# to give an optical depth of order 1 in the arbitrary demo units used
# here. When you later plug in real HITRAN parameters, you will need to
# convert S from HITRAN's wavenumber-based units into these frequency-
# based units - flagging that conversion is left for your radiative-
# transfer project.
tau_nu = line_strength * column_density * phi_voigt
transmission_nu = np.exp(-tau_nu)

plt.figure(figsize=(7, 6))

plt.subplot(2, 1, 1)
plt.plot(detuning_grid, tau_nu, color="tab:red")
plt.ylabel("Optical depth  tau_nu")
plt.title("Radiation-transfer demo: optical depth and transmission")

plt.subplot(2, 1, 2)
plt.plot(detuning_grid, transmission_nu, color="tab:green")
plt.xlabel(f"Frequency detuning  (nu - nu0)  [{freq_unit_name}]")
plt.ylabel("Transmission  exp(-tau_nu)")

plt.tight_layout()
plt.savefig("plot6_radiative_transfer_demo.png", dpi=130)
plt.close()


# ------------------------------------------------------------------------
# 17. PRINT IMPORTANT CALCULATED QUANTITIES
# ------------------------------------------------------------------------
fwhm_gaussian = 2.0 * np.sqrt(np.log(2.0)) * delta_nu_D   # Gaussian FWHM
fwhm_lorentzian = 2.0 * gamma                             # Lorentzian FWHM

# Empirical Voigt FWHM approximation (Olivero & Longbothum 1977, JQSRT 17,
# 233-236), accurate to about 0.02% over the full range of a:
fwhm_voigt_approx = 0.5346 * fwhm_lorentzian + np.sqrt(
    0.2166 * fwhm_lorentzian ** 2 + fwhm_gaussian ** 2
)

print("---------------------------------------------------------")
print("SUMMARY OF CALCULATED QUANTITIES")
print("---------------------------------------------------------")
print(f"nu0                        = {nu0:.4e} Hz")
print(f"temperature                = {temperature:.1f} K")
print(f"mass                       = {mass:.4e} kg  ({mass_amu} amu)")
print(f"delta_nu_D (Doppler width) = {delta_nu_D:.4e} Hz")
print(f"gamma (Lorentzian HWHM)    = {gamma:.4e} Hz")
print(f"a = gamma/delta_nu_D       = {a_parameter:.4f}")
print(f"Gaussian FWHM              = {fwhm_gaussian:.4e} Hz")
print(f"Lorentzian FWHM            = {fwhm_lorentzian:.4e} Hz")
print(f"Voigt FWHM (approx., O&L)  = {fwhm_voigt_approx:.4e} Hz")
print(f"Integral of Voigt profile  = {voigt_integral:.6f}  (normalization check)")
print(f"Peak optical depth tau_nu  = {tau_nu.max():.4f}")
print(f"Minimum transmission       = {transmission_nu.min():.4f}")
print("---------------------------------------------------------")
