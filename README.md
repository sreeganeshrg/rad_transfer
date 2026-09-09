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
Molecular nitrogen, **N₂**, is a homonuclear diatomic molecule. When electromagnetic radiation passes through N₂ gas, certain wavelengths are absorbed because the photon energy matches the energy difference between allowed molecular energy states. The resulting variation of absorption with wavelength is called the **N₂ absorption spectrum**.
##  Photon Absorption

Absorption of a photon promotes the molecule from an initial electronic state $i$ to an excited state $f$ according to

$$
\mathrm{N_2(X^1\Sigma_g^+)} + h\nu
\rightarrow
\mathrm{N_2^*}.
$$

The absorbed photon energy is

$$
E_\gamma = h\nu = \frac{hc}{\lambda}.
$$

Energy conservation requires

$$
E_f-E_i=\frac{hc}{\lambda_{\mathrm{abs}}}.
$$



Unlike a simple atomic spectrum, **molecular nitrogen produces bands** because a molecule has several types of energy:

$$
E=E_{\text{electronic}}+E_{\text{vibrational}}+E_{\text{rotational}}
$$

So an N₂ spectrum can contain many closely spaced lines grouped into **spectral bands**.

###  Electronic, vibrational and rotational transitions

The hierarchy is approximately:

$$
E_{\text{electronic}} \gg E_{\text{vibrational}} \gg E_{\text{rotational}}
$$

Thus:

- **Electronic transition** → produces a band system.
- **Vibrational transition** → divides an electronic band into vibrational bands.
- **Rotational transition** → produces many closely spaced rotational lines within each vibrational band.
##  Photoabsorption Cross Section

The fundamental quantity used in the simulation is the photoabsorption cross section,

$$
\sigma_{\mathrm{abs}}(\lambda)
$$

with units of cm² molecule⁻¹.

It represents the wavelength-dependent probability of interaction between an incident photon and an N₂ molecule. A large value of $\sigma_{\mathrm{abs}}$ therefore corresponds to strong absorption.

For a homogeneous gas of number density $n$, the attenuation of radiation is described by

$$
\frac{dI_\lambda}{dx} = -n\sigma_{\mathrm{abs}}(\lambda)I_\lambda
$$

Integration gives the Beer–Lambert relation,

$$
I_\lambda(L) = I_\lambda(0)\exp\left[-n\sigma_{\mathrm{abs}}(\lambda)L\right]
$$

The optical depth is therefore

$$
\tau_\lambda = nL\sigma_{\mathrm{abs}}(\lambda)
$$

and

$$
I_\lambda = I_{\lambda,0}e^{-\tau_\lambda}
$$

This relation provides the direct connection between the experimental N₂ cross-section data and the radiation-transfer simulation.
## Experimental Dataset

The experimental data used in this work are the N₂ photoabsorption cross sections reported by Souza and Srivastava (1994). The dataset corresponds to nitrogen at

$$
T = 298\ \mathrm{K}
$$

and covers the wavelength range

$$
25\ \mathrm{nm} \leq \lambda \leq 226\ \mathrm{nm}.
$$

The spectrum was obtained from an electron-energy-loss spectrum (EELS) measured using a crossed electron-beam–molecular-beam geometry. The EELS spectrum was converted into a photoabsorption spectrum, and the relative absorption cross sections were normalized to the N₂ absorption cross section at 58.6 nm.

The numerical values used in this work are obtained from the digitized data of Fig. 5 of Souza and Srivastava (1994). The resulting quantity is the wavelength-dependent molecular photoabsorption cross section,

$$
\sigma_{\mathrm{abs}}(\lambda)
\qquad
[\mathrm{cm^2\,molecule^{-1}}].
$$
## Photoabsorption Cross-Section Spectrum

##  Photoabsorption Cross-Section Spectrum

<p align="center">
  <img src="Lab%202/Initial%20phase.png" width="850">
</p>
  <b>Figure 1:</b> Photoabsorption cross section of molecular nitrogen, N₂, at 298 K over the wavelength range 25–226 nm. The spectrum shows a strong increase in absorption in the vacuum-ultraviolet region, with pronounced absorption features near 80–100 nm.
</p>

The wavelength dependence of the measured absorption cross section is shown in Fig. 1. The spectrum exhibits a strongly wavelength-dependent absorption probability, with the largest cross sections occurring in the approximately 80–100 nm region. The variation of the absorption cross section, $\sigma_{\mathrm{abs}}$, reflects the energy dependence of the electronic excitation probability of N₂.
  <b>Figure 1:</b> Photoabsorption cross section of molecular nitrogen, N₂, at 298 K over the wavelength range 25–226 nm. The spectrum shows a strong increase in absorption in the vacuum-ultraviolet region, with pronounced absorption features near 80–100 nm.
</p>

The wavelength dependence of the measured absorption cross section is shown in Fig. 1. The spectrum exhibits a strongly wavelength-dependent absorption probability, with the largest cross sections occurring in the approximately 80–100 nm region. The variation of the absorption cross section, $\sigma_{\mathrm{abs}}$, reflects the energy dependence of the electronic excitation probability of N₂.

##  Region-by-region physical interpretation

| Wavelength | Photon energy | Physical interpretation | Transition picture |
|---|---|---|---|
| 25–~75 nm | 49.6–16.5 eV | High-energy electronic absorption; continuum/photoionization becomes important as photon energy approaches and exceeds the ionization threshold. | Electronic excitation and ionization-related absorption; not a rotational/vibrational spectrum. |
| ~75–80 nm | 16.5–15.5 eV | Near the N₂ first-ionization threshold (≈15.6 eV, ≈79.5 nm). | Bound electronic absorption approaches the ionization continuum; photoionization can begin above threshold. |
| ~80–100 nm | 15.5–12.4 eV | Strongest structured absorption region in the supplied graph. | Ground-state X¹Σg⁺ → singlet ungerade electronic states, including b¹Πu, b′¹Σu⁺ and Rydberg-series states such as cₙ¹Πu / cₙ′¹Σu⁺; vibronic structure is superposed. |
| ~100–120 nm | 12.4–10.3 eV | Cross section drops strongly, but electronic/vibronic absorption remains possible. | Weak/structured electronic systems; detailed assignment requires higher-resolution data. |
| ~120–180 nm | 10.3–6.9 eV | Weak long-wavelength absorption systems. | Includes the LBH system: X¹Σg⁺ → a¹Πg in absorption; reverse a¹Πg → X¹Σg⁺ in emission. |
| ~180–226 nm | 6.9–5.5 eV | Very weak absorption in the supplied dataset. | Weak electronic/vibronic transitions; do not assign every small feature to a named transition without high-resolution data. |

The wavelength-to-energy relation used above is

$$
E(\mathrm{eV}) \approx \frac{1239.84}{\lambda(\mathrm{nm})}
$$

This makes it clear why the 25–226 nm dataset probes electronic rather than rotational or ordinary vibrational excitation.
