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

## 1.1 Photon Absorption

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

### 4. Electronic, vibrational and rotational transitions

The hierarchy is approximately:

$$
E_{\text{electronic}} \gg E_{\text{vibrational}} \gg E_{\text{rotational}}
$$

Thus:

- **Electronic transition** → produces a band system.
- **Vibrational transition** → divides an electronic band into vibrational bands.
- **Rotational transition** → produces many closely spaced rotational lines within each vibrational band.
- ## 1.2 Photoabsorption Cross Section

The fundamental quantity used in the simulation is the photoabsorption cross section,

$$
\sigma_{\mathrm{abs}}(\lambda)
$$

with units of $\mathrm{cm^2\,molecule^{-1}}$.

It represents the wavelength-dependent probability of interaction between an incident photon and an N₂ molecule. A large $\sigma_{\mathrm{abs}}$ therefore corresponds to strong absorption.

For a homogeneous gas of number density $n$, the attenuation of radiation is described by

$$
\frac{dI_\lambda}{dx}
=
-n\sigma_{\mathrm{abs}}(\lambda)I_\lambda.
$$

Integration gives the Beer–Lambert relation,

$$
I_\lambda(L)
=
I_\lambda(0)
\exp\left[-n\sigma_{\mathrm{abs}}(\lambda)L\right].
$$

The optical depth is therefore

$$
\tau_\lambda
=
nL\sigma_{\mathrm{abs}}(\lambda),
$$

and

$$
I_\lambda
=
I_{\lambda,0}e^{-\tau_\lambda}.
$$

This relation provides the direct connection between the experimental N₂ cross-section data and the radiation-transfer simulation.




