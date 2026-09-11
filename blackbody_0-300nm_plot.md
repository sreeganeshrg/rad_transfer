# Understanding the Blackbody Photon Simulation, Line by Line

This is a companion to `blackbody_photon_mc.py`. It assumes no prior coding background. Every function in the script gets its own section: what it does, why it exists, and exactly how each line of Python works. The final section maps out how to grow this code into a full atmospheric radiative transfer simulation.

---

## Part 0 — The big picture, before any code

Forget Python for a second. Here's the whole idea in plain English, using an analogy.

**Imagine a dartboard shaped like the blackbody spectrum curve** — wide and short near the tails, tall in the middle where the peak is. If you throw darts completely at random *at the board's outline* (not uniformly at the wall), you'd expect more darts to land under the tall part of the curve than under the short parts, simply because there's more "target" there. That's exactly what we want: photons should land at wavelengths in proportion to how much of the Planck spectrum lives there.

The trick the code uses to make this happen is called **inverse-transform sampling**, and it works like this:

1. Take the theoretical curve (photon-number density vs. wavelength) and turn it into a **cumulative** curve — "what fraction of all photons have a wavelength shorter than this?" This cumulative curve always starts at 0 and climbs monotonically to 1.
2. Pick a **uniform random number between 0 and 1** — every value in that range is equally likely, no shape to it at all.
3. Ask: "at what wavelength does the cumulative curve cross this random value?" That wavelength is your simulated photon.

Because the cumulative curve climbs *fast* where the original curve is *tall* (the peak) and climbs *slowly* where the original curve is short (the tails), a uniform random number is much more likely to land in the "fast-climbing" region — which automatically reproduces the correct shape, without ever needing to draw non-uniform random numbers directly. This one trick (steps 1–3) is the entire engine of the simulation. Everything else in the code is bookkeeping to make those three steps numerically accurate and fast.

The six-stage pipeline you just saw in the diagram is this idea, spelled out as code:

1. Load in the physical constants and pick a temperature.
2. Build a numerical grid of wavelengths to work on (a table of x-values).
3. Compute how much "photon density" sits at each grid point (the curve itself).
4. Turn that curve into a probability density (rescale so the area is exactly 1), then integrate it into the cumulative curve.
5. Draw uniform random numbers and read off wavelengths from the cumulative curve (steps 2–3 of the dartboard idea, done all at once for thousands of photons).
6. Check the results against known physics and draw plots.

---

## Part 1 — Setup: imports and the module docstring

```python
"""
Monte Carlo simulation of individual photons sampled from a blackbody
(Planck) spectrum.
...
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.special import zeta
```

**The triple-quoted string at the top** (`""" ... """`) is a *docstring* — a comment that describes what the file does. Python ignores it when running the code; it's there purely for humans (and tools that auto-generate documentation).

**`import numpy as np`**: NumPy is a library — someone else's pre-written code — that adds fast array math to Python. `as np` just gives it a short nickname, so instead of typing `numpy.array(...)` everywhere, you type `np.array(...)`. This is a universal convention; virtually every Python scientific script does this exact line.

**`import matplotlib.pyplot as plt`**: Matplotlib is the plotting library. `plt` is its standard nickname.

**`from scipy.optimize import brentq`**: SciPy is a library of more specialized scientific tools built on top of NumPy. This line says "from the `optimize` sub-part of SciPy, grab just the one tool called `brentq`, and let me refer to it directly as `brentq`" (rather than the longer `scipy.optimize.brentq`). `brentq` is a **root finder** — given a mathematical function, it numerically finds the input value where that function equals zero. We use it later to solve equations like "$x = 4(1-e^{-x})$" that have no clean algebraic solution.

**`from scipy.special import zeta`**: grabs the Riemann zeta function, $\zeta(s)$, which shows up in the exact formula for mean photon energy ($\zeta(3) \approx 1.20206$).


---

## Part 2 — Stage 1: physical constants

```python
h   = 6.62607015e-34    # Planck constant,        J s
c   = 2.99792458e8      # speed of light in vacuum, m/s
k_B = 1.380649e-23      # Boltzmann constant,     J/K

c2 = h * c / k_B
```

**`h   = 6.62607015e-34`**: this creates a variable named `h` and assigns it the value $6.62607015\times10^{-34}$. The `e-34` is scientific notation: `e` means "times ten to the power of," so `6.626e-34` means $6.626\times10^{-34}$. Python understands this notation natively — no special library needed for it.

Why these exact numbers, with no uncertainty stated? Since 2019, these three constants are **exact by definition** in the SI system of units — they used to be measured (with tiny error bars) and everything else was built from them; now it's the other way around. So there's no "more precise value" to look up — this literally *is* the definition.

**`c2 = h * c / k_B`**: an ordinary variable assignment where the value comes from a formula instead of a literal number. `*` is multiplication, `/` is division — standard math operators, evaluated left to right like you'd expect from a calculator. This particular combination, $c_2 = hc/k_B$, is called the "second radiation constant." It shows up so often in the Planck-law equations that it's worth giving its own name rather than writing `h * c / k_B` over and over.

## Part 3 — Stage 2: temperature

```python
T = 5000.0   # kelvin
```

Just a single number. Change this line and re-run the script, and every plot and statistic downstream updates automatically — this is the *only* place temperature is hard-coded; every function below takes temperature as an input rather than assuming this value, which is what lets the later "compare four temperatures" plot work by simply calling the same functions four times with different numbers.


---

## Part 4 — The Planck-law functions

```python
def planck_B_lambda(wavelength_m, temperature_K):
    x = c2 / (wavelength_m * temperature_K)
    return (2.0 * h * c**2 / wavelength_m**5) / np.expm1(x)


def planck_B_nu(freq_Hz, temperature_K):
    x = h * freq_Hz / (k_B * temperature_K)
    return (2.0 * h * freq_Hz**3 / c**2) / np.expm1(x)
```

**`def planck_B_lambda(wavelength_m, temperature_K):`** — `def` means "define a function." `planck_B_lambda` is the name we're giving it (names are our choice; this one is descriptive on purpose). The parentheses list its **parameters** — the inputs it expects: a wavelength in meters, and a temperature in kelvin. The colon starts the function's body, which must be indented (Python uses indentation, not curly braces, to mark "this code belongs to this function").

**`x = c2 / (wavelength_m * temperature_K)`** — inside the function, this computes the dimensionless variable $x = hc/(\lambda k_B T)$ using the `c2` shortcut from Stage 1. Note: `wavelength_m` and `temperature_K` are the *parameter names*, not global variables — whatever value gets passed in when the function is *called* is what fills these in. This is the whole point of writing a function once: the same four lines of code correctly compute $B_\lambda$ whether you call it with `wavelength_m=500e-9, temperature_K=5000` or with a million different wavelengths at once (more on that below).

**`return (2.0 * h * c**2 / wavelength_m**5) / np.expm1(x)`** — `return` hands a value back to whoever called the function. `**` is Python's "raised to the power of" operator (`c**2` means $c^2$, `wavelength_m**5` means $\lambda^5$). This line is a direct, literal translation of the Planck's-law formula:
$$B_\lambda = \frac{2hc^2}{\lambda^5}\cdot\frac{1}{e^x - 1}$$

**Why `np.expm1(x)` instead of `np.exp(x) - 1`?** This is a subtle but important numerical-stability point. When $x$ is very small (long wavelengths), $e^x$ is very close to 1 — something like 1.0000001. Computers store numbers with a fixed number of digits, so subtracting 1 from 1.0000001 the naive way can lose almost all of your precision (this is called "catastrophic cancellation"). `np.expm1(x)` computes $e^x-1$ using a method that sidesteps this — accurate to full precision even when $x$ is tiny. Concretely: for $x = 10^{-10}$, `np.exp(x) - 1` returns something noticeably wrong due to rounding, while `np.expm1(x)` returns the correct answer to 15+ decimal digits.

**One quiet but critical detail**: nothing in this function says "loop over each wavelength." If `wavelength_m` is a single number, you get a single number back. If `wavelength_m` is a NumPy array of 200,000 wavelengths, every operation (`/`, `**`, `np.expm1`) is automatically applied to *all 200,000 values at once*, and you get an array of 200,000 results back — no `for` loop anywhere. This is called **vectorization**, and it's the single biggest reason NumPy code is both short and fast. Writing an explicit loop would work too, but would be 10–100x slower and several lines longer.

`planck_B_nu` is the exact same idea in the frequency domain, using $x = h\nu/(k_BT)$ and $B_\nu = \frac{2h\nu^3}{c^2}\cdot\frac{1}{e^x-1}$.


---

## Part 5 — The dimensionless shape function

```python
def g_shape(x):
    return x**2 / np.expm1(x)
```

This is the photon-number spectrum rewritten purely in terms of $x = hc/(\lambda k_BT)$, with all the temperature and constant factors stripped away (they don't affect *where* the tolerance cutoffs should go — see Section 5–6 of the write-up for why). We use this simplified, temperature-independent version specifically to answer the question "how wide a numerical window do I need?" — described next.

## Part 6 — Choosing the numerical window: `find_x_window`

This is the most algorithmically dense function in the whole script, so we'll go slowly.

```python
def find_x_window(tolerance, x_scan_lo=1e-8, x_scan_hi=80.0, n_scan=600_000):
    x_scan = np.logspace(np.log10(x_scan_lo), np.log10(x_scan_hi), n_scan)
    g_vals = g_shape(x_scan)
    cdf = np.concatenate(
        [[0.0], np.cumsum(0.5 * (g_vals[1:] + g_vals[:-1]) * np.diff(x_scan))]
    )
    cdf /= cdf[-1]
    x_lo = np.interp(tolerance / 2.0, cdf, x_scan)
    x_hi = np.interp(1.0 - tolerance / 2.0, cdf, x_scan)
    return x_lo, x_hi
```

**The function signature**: `tolerance` is required (no default given), while `x_scan_lo=1e-8`, `x_scan_hi=80.0`, and `n_scan=600_000` are **default arguments** — if you call `find_x_window(1e-6)` without specifying the others, Python automatically uses those defaults. You *can* override them (e.g. `find_x_window(1e-6, n_scan=1_000_000)`), but usually you won't need to. (The underscore in `600_000` is just a readability separator — Python ignores it; `600_000` and `600000` are identical to the interpreter.)

**`x_scan = np.logspace(np.log10(x_scan_lo), np.log10(x_scan_hi), n_scan)`** — `np.logspace` builds an array of numbers evenly spaced *on a logarithmic scale* between two endpoints (unlike `np.linspace`, which spaces them evenly on a normal/linear scale). Concretely, `np.logspace(np.log10(1e-8), np.log10(80), 600_000)` produces 600,000 numbers starting near $10^{-8}$ and ending at 80, where each step multiplies the previous value by roughly the same *ratio* rather than adding the same *amount*. Why log-spacing here? Because the region we care about spans many orders of magnitude (from $10^{-8}$ to 80), and we need good resolution *everywhere* in that range, not just clustered near one end — log-spacing gives that automatically. (`np.log10` is just "log base 10," used here because `logspace` expects its endpoints already in log form.)

**`g_vals = g_shape(x_scan)`** — calls the function from Part 5 on the whole array at once (vectorization again), giving 600,000 values of $g(x)=x^2/(e^x-1)$.

**The trapezoidal rule, explained before the code**: to find the area under a curve numerically, one simple and reliable method is to treat the space between each pair of adjacent grid points as a thin trapezoid (flat on the bottom, sloped on top, connecting the two curve heights), compute each trapezoid's area, and add them all up. The area of one trapezoid between points $i$ and $i+1$ is:
$$\text{area}_i = \frac{g_i + g_{i+1}}{2}\times(x_{i+1}-x_i)$$
— average height times width. Do this for every adjacent pair, and you've numerically integrated the curve.

**`0.5 * (g_vals[1:] + g_vals[:-1]) * np.diff(x_scan)`** — this one line computes *all* the trapezoid areas at once. Here's each piece:
- `g_vals[1:]` uses **slicing**: `[1:]` means "give me every element starting from index 1 to the end" — i.e., everything *except the first* element. In plain terms, this is "$g_{i+1}$ for every $i$."
- `g_vals[:-1]` means "everything *except the last* element" — this is "$g_i$ for every $i$."
- Adding these two arrays together, elementwise, gives $g_i + g_{i+1}$ for every adjacent pair simultaneously.
- `np.diff(x_scan)` computes the difference between each pair of adjacent elements — exactly $x_{i+1}-x_i$ for every $i$ — in one call.
- Multiplying by `0.5` and by the diffs gives every trapezoid's area, all 599,999 of them, in a single vectorized expression.

**`np.cumsum(...)`** — "cumulative sum." Given an array like `[1, 2, 3, 4]`, `np.cumsum` returns `[1, 3, 6, 10]` — each entry is the running total of everything before it, including itself. Applied to our trapezoid areas, this gives the *running total area under the curve up to each grid point* — which is exactly the (unnormalized) cumulative distribution function (CDF).

**`np.concatenate([[0.0], np.cumsum(...)])`** — `np.concatenate` joins arrays end to end. We need this because `np.cumsum` of the trapezoid areas gives one value *per gap between points*, i.e. one fewer value than the number of grid points — but we want a CDF value *at* every grid point, starting from exactly 0 at the very first one. So we glue a single `0.0` onto the front (written as `[0.0]`, a one-item list, so `concatenate` has two arrays to join) to make the lengths match up.

**`cdf /= cdf[-1]`** — `cdf[-1]` means "the last element of the array" (negative indices count from the end in Python: `-1` is the last item, `-2` is second-to-last, and so on). `/=` divides `cdf` by that value *and stores the result back into `cdf`* (shorthand for `cdf = cdf / cdf[-1]`). Since the running total's final value is the *total* area under the curve, dividing everything by it rescales the whole CDF so it ends exactly at 1.0 — turning "raw cumulative area" into "cumulative probability."

**`np.interp(tolerance / 2.0, cdf, x_scan)`** — this is the same inverse-lookup idea from the dartboard analogy, used here for a different purpose (finding tolerance boundaries rather than sampling photons). `np.interp(target, known_x, known_y)` asks: "if I have a table of `(known_x, known_y)` pairs, and I know a target value that falls *between* two of the `known_x` entries, what `known_y` value would I expect there, by drawing a straight line between the two nearest points?" Here we flip the usual roles: we're searching *within* the `cdf` array (treating it as the "x-axis" to search) to find where it equals `tolerance/2`, and reading off the corresponding `x_scan` value (the "y-axis" being returned). In plain English: "find the $x$ value below which only `tolerance/2` of the total probability lies." The second call does the mirror image for the *upper* tail: "find the $x$ value above which only `tolerance/2` of the probability lies."

**Why split the tolerance in half?** We're trimming *both* ends of the distribution — a little bit of negligible probability from the very-long-wavelength tail, and a little bit from the very-short-wavelength tail — and we want the *total* omitted probability (both ends combined) to equal `tolerance`, so each end gets half.

The function returns `x_lo, x_hi` — Python lets a function return multiple values at once, separated by commas; whoever calls this function can unpack them into two separate variables in one line, e.g. `a, b = find_x_window(1e-6)`.


---

## Part 7 — Stage 3: building the wavelength grid

```python
TOLERANCE = 1e-6

def build_wavelength_grid(temperature_K, tolerance=TOLERANCE, n_grid=200_000):
    x_lo, x_hi = find_x_window(tolerance)
    x_grid = np.logspace(np.log10(x_lo), np.log10(x_hi), n_grid)
    wavelength_grid = c2 / (temperature_K * x_grid)
    wavelength_grid = wavelength_grid[::-1]
    return wavelength_grid, x_lo, x_hi
```

**`TOLERANCE = 1e-6`** — written in capital letters by convention (not a Python rule, just a widely-followed style) to signal "this is a constant setting, not something that changes while the program runs."

**`x_lo, x_hi = find_x_window(tolerance)`** — calls the function from Part 6 and unpacks its two return values into two variables in one line.

**`x_grid = np.logspace(np.log10(x_lo), np.log10(x_hi), n_grid)`** — same log-spacing idea as before, but now restricted to just the window we actually need (from `x_lo` to `x_hi`), with 200,000 points across that narrower range — much better resolution than scanning the full $10^{-8}$-to-$80$ range would give.

**`wavelength_grid = c2 / (temperature_K * x_grid)`** — converts each $x$ value into a wavelength using $\lambda = c_2/(Tx)$ (the definition of $x$, solved for $\lambda$). Because `x_grid` goes from small to large, and $\lambda$ is *inversely* related to $x$, this produces a `wavelength_grid` that goes from **large to small** — backwards from what we want.

**`wavelength_grid = wavelength_grid[::-1]`** — reverses the array. The slicing notation `[start:stop:step]` — here all three parts are left implicit except the step, which is `-1`, meaning "walk through the array backwards, one element at a time." `[::-1]` is a very common Python idiom for "give me this sequence reversed." After this line, `wavelength_grid` runs from shortest to longest wavelength, ascending — which every downstream function (integration, plotting) expects.

The function returns three things: the grid itself, plus `x_lo` and `x_hi` (useful later for printing diagnostic messages about the window that was chosen).

## Part 8 — Stages 4–5: the actual photon-number spectrum

```python
def photon_number_spectrum(wavelength_m, temperature_K):
    x = c2 / (wavelength_m * temperature_K)
    return (1.0 / wavelength_m**4) / np.expm1(x)


def photon_number_spectrum_nu(freq_Hz, temperature_K):
    x = h * freq_Hz / (k_B * temperature_K)
    return freq_Hz**2 / np.expm1(x)
```

`photon_number_spectrum` is $dN/d\lambda \propto \frac{1}{\lambda^4}\cdot\frac{1}{e^x-1}$ — the photon-number density derived by dividing the energy spectrum $B_\lambda$ by the energy of a single photon at that wavelength. Notice we've dropped the constant multiplying factors that appeared in the full physics derivation (like the $2c$ in front) — since this function's *only* job is to be normalized into a probability density in the next stage, any constant multiplying every point cancels out exactly when we normalize. There's no point computing it. `photon_number_spectrum_nu` is the frequency-domain twin, used later only for finding the most-probable-frequency statistic (never for sampling itself — we always sample in wavelength and derive frequency-domain quantities from those same samples).

## Part 9 — Stages 6–7: normalizing into a PDF, then integrating into a CDF

```python
def build_pdf_and_cdf(wavelength_grid, temperature_K):
    pdf_unnormalized = photon_number_spectrum(wavelength_grid, temperature_K)
    normalization = np.trapezoid(pdf_unnormalized, wavelength_grid)
    pdf = pdf_unnormalized / normalization

    increments = 0.5 * (pdf[1:] + pdf[:-1]) * np.diff(wavelength_grid)
    cdf = np.concatenate([[0.0], np.cumsum(increments)])
    cdf /= cdf[-1]
    return pdf, cdf
```

**`pdf_unnormalized = photon_number_spectrum(wavelength_grid, temperature_K)`** — evaluates the shape function on every point of our finished wavelength grid (again, vectorized — one call computes all 200,000 values).

**`normalization = np.trapezoid(pdf_unnormalized, wavelength_grid)`** — `np.trapezoid(y_values, x_values)` is NumPy's built-in trapezoidal-rule integrator — it does exactly the arithmetic we spelled out by hand in Part 6 (`0.5*(y[i]+y[i+1])*(x[i+1]-x[i])`, summed over all $i$), but as a single, well-tested library call rather than hand-written code. We use it here to compute the *total area* under our unnormalized curve.

**`pdf = pdf_unnormalized / normalization`** — dividing every point of the curve by its own total area rescales it so the *new* total area is exactly 1. This is what "normalizing" means, concretely: **a probability density is defined by having a total area of exactly 1 underneath it.** Before this line, `pdf_unnormalized` was just "relative photon counts" (correct shape, arbitrary scale); after it, `pdf` is an honest probability density you could integrate over any sub-range to get "the probability a photon falls in that range."

**The CDF construction** (`increments`, `np.concatenate`, `cdf /= cdf[-1]`) is line-for-line the same technique explained in full in Part 6 — cumulative trapezoidal integration, glue a leading zero on, then rescale so it ends at 1. We're just applying the identical recipe to the *real*, temperature-specific wavelength grid now, instead of the abstract $x$-space scan used earlier only to pick tolerance boundaries.

Why write this logic twice (once in `find_x_window`, once here) instead of one shared helper? In this case it was simpler to keep them separate since they operate on different variables ($x$ vs. $\lambda$) for different purposes (choosing boundaries vs. building the actual sampling table) — but you're right to notice the duplication; in a larger codebase you'd likely factor "cumulative trapezoidal integrate, then normalize to 1" into its own small helper function and call it from both places.


---

## Part 10 — Stages 8–9: the actual Monte Carlo sampling step

This is the heart of the whole simulation — everything before this was preparation.

```python
def sample_photon_wavelengths(n_photons, wavelength_grid, cdf, rng):
    u = rng.uniform(0.0, 1.0, n_photons)
    return np.interp(u, cdf, wavelength_grid)
```

**`rng.uniform(0.0, 1.0, n_photons)`** — `rng` is a random number generator object (created elsewhere, see Part 13). `.uniform(low, high, size)` draws `size` random numbers, each independently and equally likely to be anywhere between `low` and `high`. So `rng.uniform(0.0, 1.0, 100_000)` gives you an array of 100,000 numbers, each uniformly random in $[0,1]$ — step 2 of the dartboard analogy from Part 0, done all at once.

**`np.interp(u, cdf, wavelength_grid)`** — this is the payoff line, and it's worth working through with real numbers. Suppose (toy example, not the real physics) our CDF table looked like this:

| wavelength (nm) | cumulative probability |
|---|---|
| 500 | 0.10 |
| 600 | 0.40 |
| 700 | 0.75 |
| 800 | 0.95 |

If one of our random draws happened to be `u = 0.60`, `np.interp` looks at the `cdf` column, finds that 0.60 falls between the table rows for 600 nm (cdf 0.40) and 700 nm (cdf 0.75), and does a straight-line (linear) interpolation between those two points to estimate the wavelength: since 0.60 is $\frac{0.60-0.40}{0.75-0.40} \approx 0.57$ of the way from 0.40 to 0.75, the estimated wavelength is $600 + 0.57\times(700-600) \approx 657$ nm. That number — 657 nm in this toy example — is one simulated photon.

Do this for 100,000 different random `u` values simultaneously (which `np.interp` happily does when you hand it a whole array instead of one number), and you get 100,000 simulated photon wavelengths in one line of code, each one automatically respecting the correct Planck-spectrum proportions — *because* the table (`cdf`, `wavelength_grid`) it's interpolating against was built directly from the physics in Parts 5–9.


---

## Part 11 — Stage 12: statistics with `summarize`

```python
def summarize(samples, wavelength_grid, pdf, temperature_K, label=""):
    photon_energies = h * c / samples
    mean_energy = photon_energies.mean()
    mean_energy_theory = (np.pi**4 / (30.0 * zeta(3))) * k_B * temperature_K

    mean_wavelength = samples.mean()
    i_peak = np.argmax(pdf)
    mode_wavelength = wavelength_grid[i_peak]

    x_peak_photon_nu = brentq(lambda x: 2 * (1 - np.exp(-x)) - x, 0.3, 5.0)
    mode_frequency = x_peak_photon_nu * k_B * temperature_K / h
    ...
    return {"mean_energy": mean_energy, ...}
```

**`photon_energies = h * c / samples`** — `samples` is the whole array of simulated wavelengths from Part 10. Dividing `h*c` by the *entire array* at once (vectorization, again) applies $E = hc/\lambda$ to every simulated photon simultaneously, giving an array of individual photon energies.

**`.mean()`** — a method (a function attached to an object, called with a dot) built into every NumPy array. `photon_energies.mean()` adds up all the energies and divides by how many there are — the ordinary arithmetic mean — in one call, instead of writing a manual running-total loop.

**`mean_energy_theory = (np.pi**4 / (30.0 * zeta(3))) * k_B * temperature_K`** — a direct transcription of the closed-form result $\langle E\rangle = \frac{\pi^4}{30\zeta(3)}k_BT$. `np.pi` is just NumPy's stored value of $\pi$; `zeta(3)` calls the imported SciPy function to compute $\zeta(3)\approx1.2021$. This line exists purely so we can compare our *simulated* mean energy against the *known correct answer* — the validation step.

**`i_peak = np.argmax(pdf)`** — `np.argmax` returns the *index* (position) of the largest value in an array, not the value itself. If `pdf` were `[0.1, 0.4, 0.9, 0.3]`, `np.argmax(pdf)` would return `2` (the position of `0.9`, counting from 0).

**`mode_wavelength = wavelength_grid[i_peak]`** — array indexing: `wavelength_grid[2]` means "the 3rd element of `wavelength_grid`" (Python counts array positions starting from 0, not 1 — so index `2` is actually the *third* item). Combined with the line above, these two lines together answer "at which wavelength does the probability density reach its highest point?" — i.e., the most probable wavelength.

**`brentq(lambda x: 2 * (1 - np.exp(-x)) - x, 0.3, 5.0)`** — this is root-finding, mentioned back in Part 1. We want to solve the equation $x = 2(1-e^{-x})$ — a "transcendental" equation with no algebraic formula for $x$. `brentq` finds it numerically instead: you give it a function that equals zero exactly where your equation is satisfied (here, rearranged to `2*(1-exp(-x)) - x`), plus a bracket `(0.3, 5.0)` you're confident the answer sits inside, and it searches that bracket until it converges on the root.

**`lambda x: 2 * (1 - np.exp(-x)) - x`** is a **lambda function** — a tiny, throwaway, unnamed function defined inline, for cases where writing a full separate `def` block would be overkill for something used once. It means exactly the same thing as:
```python
def temp_function(x):
    return 2 * (1 - np.exp(-x)) - x
```
just written on one line, without needing to invent and remember a name for it.

**The dictionary return value**: `{"mean_energy": mean_energy, ...}` is a Python **dictionary** — a collection of `key: value` pairs, where you look things up by name (`key`) rather than by numeric position. Whoever calls `summarize()` can then do `result["mean_energy"]` to get that specific number back out, which is more readable than remembering "the third thing this function returns."

**f-strings** (used throughout the `print(...)` lines, e.g. `f"mean energy: {mean_energy:.4e} J"`): the `f` before the quote marks means "this string can have `{...}` placeholders that get filled in with actual variable values." The `:.4e` part is a *format spec* — it means "display this number in scientific notation with 4 digits after the decimal point." This is purely for making printed output readable; it has no effect on the underlying numbers.

## Part 12 — `fraction_in_range`: boolean masks

```python
def fraction_in_range(samples, lam_min, lam_max):
    inside = (samples >= lam_min) & (samples <= lam_max)
    return inside.mean()
```

**`samples >= lam_min`** — comparing an entire array to a single number produces a new array of `True`/`False` values, one per element — `True` wherever that element satisfies the comparison, `False` where it doesn't. This is called a **boolean mask**.

**`&`** — "and," applied elementwise between two boolean arrays. `(samples >= lam_min) & (samples <= lam_max)` gives `True` only where *both* conditions hold — i.e., wherever a sampled wavelength falls inside `[lam_min, lam_max]`.

**`inside.mean()`** — here's a small trick worth calling out explicitly: in Python, `True` behaves like the number 1 and `False` behaves like 0. So averaging an array of `True`/`False` values is exactly the same as computing "what fraction of these are `True`" — which is exactly "what fraction of simulated photons fell in this wavelength range." No explicit counting loop needed.


---

## Part 13 — The plotting functions and matplotlib basics

All four plotting functions follow the same pattern, so understanding one means understanding all four:

```python
def plot_histogram_vs_theory(samples, wavelength_grid, pdf, temperature_K, outpath, n_bins=120):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(samples * 1e9, bins=n_bins, density=True, ...)
    ax.plot(wavelength_grid * 1e9, pdf * 1e-9, ...)
    ax.set_xlabel("wavelength  (nm)")
    ax.set_ylabel("probability density")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
```

**`fig, ax = plt.subplots(figsize=(7, 4.5))`** — matplotlib's standard opening move. `fig` (the "figure") is the whole image/canvas; `ax` (the "axes") is the actual plotting area with its x- and y-axis, where data actually gets drawn. `figsize=(7, 4.5)` sets the image dimensions in inches (width, height). Almost every matplotlib script starts this way.

**`ax.hist(samples * 1e9, bins=n_bins, density=True, ...)`** — draws a histogram. `samples * 1e9` converts our wavelengths from meters to nanometers just for display (multiplying an entire array by a single number, again vectorized — every element gets multiplied). `bins=n_bins` sets how many bars to divide the range into. `density=True` is important: it rescales the histogram so its *own* total area is 1 — matching the scale of our theoretical probability density curve, so the two can be overlaid meaningfully on the same axes (without this, the histogram would just show raw counts, whose height depends on how many photons we simulated and how wide we made the bins — not directly comparable to a probability density).

**`ax.plot(wavelength_grid * 1e9, pdf * 1e-9, ...)`** — draws the theoretical curve as a line on the *same* axes as the histogram. `pdf * 1e-9` here isn't a units bug — it's a unit *conversion*: our `pdf` was built as "probability per meter" (since wavelength was in meters throughout the physics), but we're displaying the x-axis in nanometers, so the y-axis needs converting too, from "per meter" to "per nanometer," which is a factor of $10^{-9}$.

**`ax.set_xlabel(...)`, `ax.set_ylabel(...)`, `ax.legend()`** — labeling. `.legend()` reads the `label=...` argument you passed to `.hist()` and `.plot()` earlier and draws a small key in the corner of the plot.

**`fig.tight_layout()`** — automatically adjusts spacing so labels and titles don't get clipped or overlap the edges of the image.

**`fig.savefig(outpath, dpi=150)`** — writes the figure to disk as an actual image file (a PNG, in our case) at the given path, at 150 dots per inch (a resolution setting — higher means a sharper, larger file).

**`plt.close(fig)`** — releases the figure from memory. This matters specifically because our script makes *several* plots in a row (four different functions, one temperature-comparison loop) — without explicitly closing each figure, matplotlib keeps every single one open in memory simultaneously, which can quietly use up a lot of RAM if you're generating many plots in one script run.

The other three plotting functions (`plot_theoretical_distribution`, `plot_N_convergence`, `plot_temperature_comparison`) use these exact same building blocks — `subplots`, `.hist`/`.plot`, labels, `savefig`, `close` — just with different data fed in, and `plot_N_convergence` additionally uses `plt.subplots(1, len(N_list), ...)` to create *several* side-by-side axes in one figure (a grid of small plots) rather than just one.

---

## Part 14 — The `if __name__ == "__main__":` block, and execution order

```python
if __name__ == "__main__":
    rng = np.random.default_rng(42)
    wavelength_grid, x_lo, x_hi = build_wavelength_grid(T)
    pdf, cdf = build_pdf_and_cdf(wavelength_grid, T)
    samples_main = sample_photon_wavelengths(100_000, wavelength_grid, cdf, rng)
    stats_main = summarize(samples_main, wavelength_grid, pdf, T, label="main run")
    ...
```

**`if __name__ == "__main__":`** — this looks strange the first time you see it, so here's exactly what it means. Every Python file has a hidden built-in variable called `__name__`. If you *run* this file directly (`python3 blackbody_photon_mc.py`), Python automatically sets `__name__` equal to the string `"__main__"` for that file — so this `if` is `True`, and everything indented underneath it runs. But if you instead *import* this file from a different script (`import blackbody_photon_mc`) in order to reuse its functions, `__name__` is set to the file's actual module name instead, so this block is skipped entirely — you get access to all the function *definitions* above it, without the script automatically re-running its whole simulation and popping up four plots the moment you import it. This is standard practice in essentially every serious Python script.

**`rng = np.random.default_rng(42)`** — creates the random number generator object used throughout the script. `42` is the **seed** — a starting point for the pseudo-random number sequence. Passing the same seed always produces the *exact same* sequence of "random" numbers every time you run the script — genuinely useful for debugging and for letting someone else reproduce your exact results, even though the numbers still look statistically random. If you wanted a *different* random outcome each run, you'd either use a different seed or omit it (`np.random.default_rng()` picks an unpredictable seed on its own).

**The execution order** below this line is just a straight-line script: build the grid, build the PDF/CDF, draw the main batch of samples, summarize them, run the Wien's-law check, generate all four plots, then loop over four temperatures repeating the sampling-and-summarizing step for each. Nothing here is more complicated conceptually than "call the functions we already defined, in the order we need their outputs" — all the real logic lives in the function definitions above it.


---

## Part 15 — From here to atmospheric radiative transfer

This is the actual reason for building the blackbody simulator first: it is the source term for the bigger problem. Here's exactly what carries over unchanged, what has to change, and a concrete sequence of small steps between here and there.

### What you already have, and can reuse almost as-is

- **The blackbody sampler itself.** The Sun (roughly a 5778 K blackbody) or the Earth's surface (roughly 288 K, emitting thermal infrared) are both blackbody-ish sources. `build_wavelength_grid`, `build_pdf_and_cdf`, and `sample_photon_wavelengths` are exactly the machinery you'll use to generate a physically correct population of source photons — you'll literally call these same three functions, just with a different `T`, as the very first step of an atmospheric simulation ("here are 100,000 photons leaving the top of the atmosphere or the surface, with wavelengths sampled correctly for that source's temperature").
- **The inverse-transform-sampling *technique*.** This turns out to be exactly the tool you need for the new physics below too — not a new concept, the same one applied to a different distribution.
- **The general pattern**: normalize a distribution, build a CDF, invert it with `np.interp`. You'll reuse this pattern repeatedly.

### What's new: the atmosphere adds a "did this photon survive its journey?" question

The blackbody simulation asks one question per photon: *what wavelength is it?* Atmospheric radiative transfer asks a second, sequential question for every photon, once it has a wavelength: *as it travels through the atmosphere, does it get absorbed, scattered, or does it pass straight through?* This turns each photon's simulation from "one random draw" into "one random draw for wavelength, then a chain of random draws as it moves through the atmosphere," so the code's shape will change from "generate an array of samples all at once" to "track each photon step by step until something happens to it."

**Beer–Lambert absorption.** The probability a photon survives passing through a layer of gas without being absorbed is
$$P(\text{survive}) = e^{-\tau}$$
where $\tau$ (optical depth) is the absorption coefficient integrated along the photon's path through that layer. Here's the genuinely nice surprise: the standard Monte Carlo way to decide *how far* a photon travels before it interacts is to draw a uniform random number $u$ and set
$$\tau_{\text{travelled}} = -\ln(u)$$
This is inverse-transform sampling again — applied to the exponential distribution instead of the Planck distribution, but it's the exact same idea as `np.interp(u, cdf, wavelength_grid)`: a uniform random number, mapped through an inverted cumulative distribution, to get a physically meaningful quantity (there, a wavelength; here, a distance in optical-depth units). You already know how to do this.

**Wavelength-dependent absorption.** Real gases (water vapor, CO₂, ozone, and so on) absorb much more strongly at some wavelengths than others (molecular absorption lines/bands). Since your photon already has a wavelength assigned at emission (from the blackbody sampler), you'll look up an absorption cross-section $\sigma(\lambda)$ for that specific wavelength, multiply by the local number density of the absorbing gas to get an absorption coefficient, and that's what feeds into $\tau$. This is genuinely new data you'll need from outside your own code — high-resolution absorption line data (the HITRAN database is the standard source) or simplified band models to start with.

**Layered atmosphere structure.** Rather than "empty space," a photon now moves through a stack of horizontal layers (a "plane-parallel atmosphere," the standard simplifying assumption), each with its own temperature, pressure, and gas composition (density generally falls off roughly exponentially with altitude). Each layer contributes its own piece to the total optical depth along the photon's path. This is naturally represented as an array too — e.g. `layer_altitudes`, `layer_density`, `layer_temperature` — so a lot of the "build an array, evaluate a function on it" style from this script carries over directly.

**Scattering.** Rayleigh scattering (by air molecules, wavelength-dependent, stronger at short wavelengths — this is why the sky is blue) and Mie scattering (by aerosols/cloud droplets, more complex angular dependence) redirect a photon's *direction* rather than destroying it. Once you add scattering, each photon needs a direction (an angle, or a 3D unit vector) as part of its state, not just a wavelength — a genuinely new piece of bookkeeping. A scattering event means: pick a new direction (sampled from that scattering process's angular probability distribution — inverse-transform sampling, again, on yet another distribution), and continue the photon's journey from there.

**Photon bookkeeping.** Concretely, each simulated photon will need to carry around: its wavelength (fixed at emission, for the simple case of pure scattering/absorption with no wavelength-shifting), its current position/altitude, and its current direction. Its journey ends when one of three things happens: it's absorbed (removed from the simulation), it exits the top of the atmosphere (counted as "escaped to space"), or it reaches the surface (counted as "transmitted" or "reflected," depending on surface properties).

### A concrete, incremental path (matching how you built this one)

1. **Single layer, absorption only, one fixed wavelength.** Just apply $e^{-\tau}$ as a survival probability for a beam of monochromatic photons through one slab of gas. This is barely more than a new small function on top of what you already have.
2. **Multiple layers, still absorption only, but wavelength sampled from your blackbody sampler.** Now each photon carries a wavelength (reuse Parts 7–10 unchanged), and you look up a (possibly simplified, e.g. two or three bands) absorption coefficient per layer per wavelength, accumulating $\tau$ layer by layer.
3. **Add Rayleigh scattering.** Now photons need a direction, and a scattering event redirects them rather than ending their journey — this is where the "one random draw and done" architecture of the current script genuinely needs to become a loop ("keep advancing this photon, layer by layer, until it's absorbed or exits").
4. **Realistic absorption data and a proper atmospheric profile** (e.g., the U.S. Standard Atmosphere for temperature/pressure/density vs. altitude, and real molecular absorption cross-sections), replacing the simplified stand-ins from steps 1–3.

Each of those is a genuinely manageable next project on its own — and at every step, the core trick (turn a physical distribution into a CDF, invert it with a uniform random draw) is one you've already built and tested here.

