# rad_transfer

This repo is for our Advanced Computational Lab 
We are simulating a black body radiation transferring through a uniform medium 

## Initial Phase

<p align="center">
  <img src="Initial%20phase.png" width="850">
</p>


## week-1##
Gr leader: Shreeganesh
contributions:

## week-2##
Gr leader: Usha
Contributions: (Sarfaraz & Asha: Theory of N2-N2 absorption and experimental data collection of all absorptions in 25nm-200nm- UV region)
Future work: implementation of these absorption spectra of N2-N2 molecules in the black-body radiation spectrum in the UV region, 
Explanation of work done:


## Basic Theory of N₂ Absorption Spectrum
# Photoabsorption Spectrum of Molecular Nitrogen

## 1. Photoabsorption of Molecular Nitrogen

Molecular nitrogen, \(N_2\), is a homonuclear diatomic molecule whose ground electronic state is

\[
X^1\Sigma_g^+
\]

The photoabsorption spectrum considered in this work corresponds to the experimentally derived \(N_2\) photoabsorption cross section at \(298\,\mathrm{K}\) over approximately \(25-226\,\mathrm{nm}\).

The experimental dataset used in this work is from Souza and Srivastava (1994). The relative absorption cross sections were obtained from an electron-energy-loss spectrum and normalized to an absolute \(N_2\) absorption cross section at \(58.6\,\mathrm{nm}\).

---

## 1.1 Photon Absorption

Absorption of a photon promotes the molecule from an initial electronic state \(i\) to an excited state \(f\):

\[
N_2(X^1\Sigma_g^+) + h\nu
\rightarrow
N_2^*
\]

The absorbed photon energy is

\[
E_\gamma = h\nu = \frac{hc}{\lambda}
\]

where:

- \(h\) is Planck's constant,
- \(\nu\) is the photon frequency,
- \(c\) is the speed of light,
- \(\lambda\) is the photon wavelength.

Energy conservation requires

\[
E_f-E_i
=
\frac{hc}{\lambda_{\mathrm{abs}}}
\]

Therefore, the wavelength corresponding to an absorption transition is

\[
\boxed{
\lambda_{\mathrm{abs}}
=
\frac{hc}{E_f-E_i}
}
\]

For practical calculations, the photon energy in electron volts can be written as

\[
\boxed{
E_\gamma(\mathrm{eV})
=
\frac{1239.84}{\lambda(\mathrm{nm})}
}
\]

Thus, the wavelength range \(25-226\,\mathrm{nm}\) corresponds approximately to

\[
5.49\,\mathrm{eV}
\leq
E_\gamma
\leq
49.59\,\mathrm{eV}
\]

---

## 1.2 Photoabsorption Cross Section

The fundamental quantity used in the simulation is the photoabsorption cross section,

\[
\boxed{
\sigma_{\mathrm{abs}}(\lambda)
}
\]

with units of

\[
\mathrm{cm^2\,molecule^{-1}}
\]

The photoabsorption cross section represents the wavelength-dependent probability or strength of interaction between an incident photon and an \(N_2\) molecule.

A large value of \(\sigma_{\mathrm{abs}}\) corresponds to strong absorption, while a small value corresponds to weak absorption.

For a homogeneous gas with number density \(n\), the attenuation of radiation is described by

\[
\frac{dI_\lambda}{dx}
=
-n\sigma_{\mathrm{abs}}(\lambda)I_\lambda
\]

where:

- \(I_\lambda\) is the spectral intensity,
- \(n\) is the \(N_2\) molecular number density,
- \(\sigma_{\mathrm{abs}}(\lambda)\) is the absorption cross section,
- \(x\) is the propagation distance.

Integration over a path length \(L\) gives the Beer--Lambert relation:

\[
\boxed{
I_\lambda(L)
=
I_\lambda(0)
\exp[-n\sigma_{\mathrm{abs}}(\lambda)L]
}
\]

The optical depth is therefore

\[
\boxed{
\tau_\lambda
=
nL\sigma_{\mathrm{abs}}(\lambda)
}
\]

and the transmitted intensity can be written as

\[
\boxed{
I_\lambda
=
I_{\lambda,0}e^{-\tau_\lambda}
}
\]

The transmission is

\[
\boxed{
T(\lambda)
=
\frac{I_\lambda}{I_{\lambda,0}}
=
e^{-\tau_\lambda}
}
\]

and the absorbed fraction is

\[
\boxed{
A(\lambda)
=
1-e^{-\tau_\lambda}
}
\]

This provides the direct connection between the experimental \(N_2\) cross-section data and the radiation-transfer simulation.

---

## 1.3 Experimental N₂ Dataset

The experimental data used in this work are the \(N_2\) photoabsorption cross sections reported by Souza and Srivastava (1994).

The dataset corresponds to nitrogen at

\[
T=298\,\mathrm{K}
\]

and covers the wavelength range

\[
25\,\mathrm{nm}
\leq
\lambda
\leq
226\,\mathrm{nm}
\]

The spectrum was obtained from an electron-energy-loss spectrum (EELS) measured using a crossed electron-beam--molecular-beam geometry.

The EELS spectrum was converted into a photoabsorption spectrum, and the relative absorption cross sections were normalized to the \(N_2\) absorption cross section at \(58.6\,\mathrm{nm}\).

The numerical data are represented as the wavelength-dependent molecular photoabsorption cross section,

\[
\sigma_{\mathrm{abs}}(\lambda)
\quad
[\mathrm{cm^2\,molecule^{-1}}]
\]

---

## 1.4 Origin of Absorption Peaks

An absorption peak occurs when the photon energy corresponds to a molecular transition having appreciable transition probability.

The electric-dipole transition strength is determined by the transition dipole moment,

\[
\boxed{
M_{if}
=
\langle\Psi_f|\mu|\Psi_i\rangle
}
\]

where:

- \(\Psi_i\) is the initial molecular state,
- \(\Psi_f\) is the final molecular state,
- \(\mu\) is the electric dipole moment operator.

A transition is electric-dipole allowed when

\[
\boxed{
M_{if}\neq0
}
\]

For a homonuclear molecule such as \(N_2\), the principal electric-dipole selection rule is

\[
\boxed{
g\leftrightarrow u
}
\]

The spin selection rule is

\[
\boxed{
\Delta S=0
}
\]

Since the ground state of \(N_2\) is

\[
X^1\Sigma_g^+
\]

strong one-photon absorption is expected for excited singlet ungerade states, particularly states having

\[
^1\Sigma_u^+
\]

and

\[
^1\Pi_u
\]

symmetry.

---

## 1.5 Why the 80--100 nm Region is Strong

The most pronounced absorption in the experimental spectrum occurs approximately in the

\[
80-100\,\mathrm{nm}
\]

region.

The corresponding photon energies are approximately

\[
12.4\,\mathrm{eV}
\lesssim
E_\gamma
\lesssim
15.5\,\mathrm{eV}
\]

This region contains strong electric-dipole-allowed transitions from the ground state

\[
X^1\Sigma_g^+
\]

to singlet ungerade excited states of predominantly

\[
^1\Sigma_u^+
\]

and

\[
^1\Pi_u
\]

character.

The excited states possess both valence and Rydberg character. Strong interactions between states having the same symmetry can redistribute their energies and oscillator strengths.

Therefore, the spectrum contains several strong and irregular features rather than a single isolated absorption line.

The strong absorption in this region results from the combination of:

\[
\boxed{
\text{Allowed symmetry}
+
\text{Large transition strength}
+
\text{Rydberg--valence interaction}
}
\]

---

## 1.6 Vibronic Structure

Each electronic state of a diatomic molecule contains vibrational and rotational levels.

Therefore, an electronic transition can be written more accurately as

\[
(e'',v'',J'')
\rightarrow
(e',v',J')
\]

where:

- \(e\) represents the electronic state,
- \(v\) is the vibrational quantum number,
- \(J\) is the rotational quantum number.

The transition energy can be approximately expressed as

\[
\boxed{
\Delta E
=
\Delta E_{\mathrm{elec}}
+
\Delta E_{\mathrm{vib}}
+
\Delta E_{\mathrm{rot}}
}
\]

Thus, one electronic transition can produce many closely spaced vibronic and rotational transitions.

When these individual transitions are not resolved experimentally, they appear as broader absorption features in the measured cross section.

Therefore, the observed spectrum is better described as an

\[
\boxed{
\text{electronic/vibronic absorption spectrum}
}
\]

rather than a purely electronic spectrum.

---

## 1.7 Rydberg States and Peak Intensities

At high excitation energies, \(N_2\) exhibits Rydberg states in which an electron occupies a high-lying orbital.

These states occur in series that converge toward ionic limits:

\[
\boxed{
E_n\rightarrow E_{\mathrm{ion}}
}
\]

Rydberg states having the same symmetry as valence states can interact with each other.

This Rydberg--valence interaction changes both the energies and the distribution of oscillator strength among the observed molecular levels.

Therefore, experimentally observed absorption peaks do not necessarily correspond to isolated transitions.

Instead, the measured cross section reflects the coupled molecular states and their transition strengths.

---

## 1.8 Ionization Threshold

At sufficiently high photon energies, the \(N_2\) molecule can undergo photoionization:

\[
\boxed{
N_2+h\nu
\rightarrow
N_2^+ + e^-
}
\]

The first ionization energy of \(N_2\) is approximately

\[
E_{\mathrm{ion}}
\approx
15.6\,\mathrm{eV}
\]

The corresponding wavelength is

\[
\lambda_{\mathrm{ion}}
=
\frac{1239.84}{15.6}
\]

giving approximately

\[
\boxed{
\lambda_{\mathrm{ion}}
\approx79.5\,\mathrm{nm}
}
\]

Thus, approximately \(79.5\,\mathrm{nm}\) provides an important physical boundary in the spectrum.

For wavelengths shorter than approximately \(79.5\,\mathrm{nm}\), the photon has sufficient energy to ionize \(N_2\).

The short-wavelength region can therefore contain contributions from both discrete electronic excitation and continuum/photoionization processes.

---

## 1.9 Physical Interpretation of the Spectrum

The wavelength dependence of the measured cross section can be interpreted as

\[
\boxed{
\sigma_{\mathrm{abs}}(\lambda)
\rightarrow
\text{molecular transition probability}
}
\]

Strong values of \(\sigma_{\mathrm{abs}}\) occur where the photon energy matches transitions with substantial oscillator strength.

Weak values can occur when:

- the available transitions have small oscillator strength,
- the transition is symmetry forbidden,
- the transition is spin forbidden,
- the photon energy lies away from strong resonances.

Therefore, the broad structure observed in the experimental spectrum results from the combined effects of

\[
\boxed{
\text{Electronic transitions}
+
\text{Vibrational structure}
+
\text{Rotational structure}
+
\text{Rydberg--valence interaction}
+
\text{Ionization/continuum effects}
}
\]

---

## 1.10 Connection to Radiation-Transfer Simulation

For radiation-transfer calculations, the experimentally measured

\[
\sigma_{\mathrm{abs}}(\lambda)
\]

is used directly in the optical-depth equation:

\[
\boxed{
\tau_\lambda
=
nL\sigma_{\mathrm{abs}}(\lambda)
}
\]

The transmitted radiation is then

\[
\boxed{
I_\lambda
=
I_{\lambda,0}e^{-\tau_\lambda}
}
\]

Therefore, the experimental molecular cross section determines the wavelength-dependent attenuation of the incident radiation.

---

## 1.11 Connection with Blackbody Radiation

In this project, blackbody radiation is used as the incident radiation source.

The photon energy is

\[
E_\gamma
=
\frac{hc}{\lambda}
\]

The blackbody photon-number distribution is proportional to

\[
S_\gamma(\lambda)
\propto
\frac{1}
{\lambda^4
\left[
\exp\left(
\frac{hc}{\lambda k_BT}
\right)-1
\right]}
\]

Photon wavelengths are sampled using a Monte Carlo method.

For each sampled photon wavelength \(\lambda_i\), the experimental \(N_2\) absorption cross section is obtained:

\[
\sigma_i
=
\sigma_{\mathrm{abs}}(\lambda_i)
\]

The optical depth for each photon is

\[
\tau_i
=
nL\sigma_i
\]

The transmission probability is

\[
\boxed{
T_i=e^{-\tau_i}
}
\]

and the absorption probability is

\[
\boxed{
P_{\mathrm{abs},i}
=
1-e^{-\tau_i}
}
\]

Thus, photons at wavelengths where \(N_2\) has a large absorption cross section have a higher probability of being absorbed.

---

## 1.12 Complete Simulation Workflow

The complete simulation can be summarized as:

\[
\boxed{
\text{Blackbody Radiation}
\rightarrow
\text{Photon-Number Distribution}
\rightarrow
\text{Monte Carlo Photon Sampling}
\rightarrow
\text{Experimental } \sigma_{\mathrm{abs}}(\lambda)
\rightarrow
\text{Optical Depth}
\rightarrow
\text{Transmission}
\rightarrow
\text{Absorbed Spectrum}
}
\]

For the incident blackbody spectrum,

\[
S_{\mathrm{BB}}(\lambda)
\]

the transmitted spectrum is

\[
\boxed{
S_{\mathrm{trans}}(\lambda)
=
S_{\mathrm{BB}}(\lambda)
\exp[-nL\sigma_{\mathrm{abs}}(\lambda)]
}
\]

The absorbed spectrum is

\[
\boxed{
S_{\mathrm{abs}}(\lambda)
=
S_{\mathrm{BB}}(\lambda)
\left[
1-\exp(-nL\sigma_{\mathrm{abs}}(\lambda))
\right]
}
\]

---

## 1.13 Summary

The \(N_2\) photoabsorption spectrum is governed by transitions from the ground electronic state

\[
X^1\Sigma_g^+
\]

to electronically excited states.

The photon energy determines the accessible molecular states, while the transition dipole moment and molecular selection rules determine their absorption strength.

The strong absorption in the approximately \(80-100\,\mathrm{nm}\) region is associated with electric-dipole-allowed transitions to singlet ungerade states of

\[
^1\Sigma_u^+
\]

and

\[
^1\Pi_u
\]

character, with significant Rydberg--valence interaction.

The experimentally measured wavelength-dependent cross section

\[
\sigma_{\mathrm{abs}}(\lambda)
\]

provides the microscopic absorption property required for the radiation-transfer simulation.

The central relationship used in the simulation is

\[
\boxed{
\tau_\lambda=nL\sigma_{\mathrm{abs}}(\lambda)
}
\]

which connects the molecular absorption data to the attenuation of blackbody radiation through the Beer--Lambert law.
