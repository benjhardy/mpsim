# Multipath SNR Oscillation and Attenuation Model

This document describes the code used in this repository to capture multipath SNR oscillations and attenuation when satellites transmit from low elevation. It is intended for porting this model to other simulations (e.g., Remcom Wireless InSite) for verification against measured GNSS data.

---

## Overview

The model treats GNSS SNR as the interference between a **direct** (line-of-sight) signal and a **reflected** (multipath) signal. At low elevation angles, the path-length difference between direct and reflected rays changes slowly with elevation, producing characteristic **oscillations** in SNR. The **attenuation** of the reflected signal is governed by Fresnel coefficients, surface roughness, and antenna gain patterns.

---

## 1. Core Physical Model

### Composite Signal

The received voltage is the phasor sum of direct and reflected signals:

```
phasor_composite = phasor_direct + phasor_reflected
```

The observed power (and thus SNR) is:

```
power_composite = |phasor_composite|² = |direct + reflected|²
```

This produces constructive and destructive interference as the relative phase changes with satellite elevation.

### Key Files

- **`lib/snr/fwd/snr_fwd.m`** – Main forward simulator entry point
- **`lib/snr/fwd/snr_fwd_combine.m`** – Combines direct and reflected phasors
- **`lib/snr/fwd/snr_fwd_direct_and_reflected.m`** – Computes both signal components

---

## 2. Interferometric Phase (Oscillation Source)

### Path-Length Difference

For a horizontal reflecting surface with antenna height `h` above the surface and satellite elevation `elev` (degrees, 0° = horizon, 90° = zenith):

```
delay = 2 * h * sin(elev)
```

Where `delay` is the extra path length (meters) of the reflected ray compared to the direct ray.

**File:** `lib/snr/util/get_reflection_delay.m`

### Phase to Phasor

The reflection introduces a phase shift:

```
phase_deg = (delay / wavelength) * 360
phasor_delay = exp(j * phase_deg * π/180)
```

The reflected phasor is multiplied by this delay phasor. As elevation changes, `delay` changes, so the phase rotates and the interference oscillates.

**File:** `lib/snr/fwd/snr_fwd_delay.m`

---

## 3. Oscillation Frequency vs. Elevation

### Why sin(elevation)?

The interferometric phase varies with elevation as:

```
φ ∝ (2h/λ) * sin(elev)
```

So the **natural independent variable** for the oscillation is **sin(elevation)**, not elevation itself. At low elevation, sin(elev) is small and changes slowly, so fringes are widely spaced and easier to resolve.

### Multipath LSSA (Least-Squares Spectral Analysis)

The code fits a sinusoid to SNR (or detrended SNR) as a function of sin(elevation):

```
freq = 2 * height / wavelength     # cycles per unit sin(elev)
period = 1 / freq                  # period in sin(elev) space
sine = sin(elev)                   # independent variable
```

**File:** `lib/snr/util/mplsqfourier/mplsqfourier.m`

The dominant oscillation has:
- **Frequency:** `f = 2h/λ` (in cycles per unit sin(elev))
- **Period in sin(elev):** `λ/(2h)`
- **Amplitude:** Depends on direct/reflected power ratio
- **Phase:** Depends on reflector height and bias

### Fringe Fitting

**File:** `lib/snr/bias/snr_bias_fit_fringe.m`

- Detrends SNR (removes polynomial trend)
- Fits sinusoid via `mplsqfourier` to extract amplitude, phase, height
- Default elevation range for fringe fitting: **`elev_lim_fringe = [0 35]°`** (low elevation)

---

## 4. Extracting the Oscillation (Multipath Modulation)

### get_multipath_modulation

**File:** `lib/snr/util/get_multipath_modulation.m`

Extracts the oscillating component from composite power:

```
composite_pow = |phasor_direct + phasor_reflected|²
trend = |direct|² + |reflected|²                    # non-oscillating part
var = composite_pow - trend                        # oscillation
den = 2 * √(direct_pow * reflected_pow)             # normalization
var_normalized = var / den                           # cosine-like [-1, 1]
```

The normalized `var` approximates `cos(φ)` where φ is the interferometric phase.

### Fringe Visibility

**File:** `lib/snr/util/get_fringe_vis.m`

```
power_interf = |phasor_reflected / phasor_direct|²
fringe_vis = 2√(power_interf) / (1 + power_interf)
```

Fringe visibility quantifies how strong the oscillation is (0 = no oscillation, 1 = full contrast).

### Fringe Deviation (Peak-to-Trough in dB)

**File:** `lib/snr/util/get_fringe_dev.m`

Computes the dB difference between constructive and destructive interference (max vs. min of the oscillation).

---

## 5. Attenuation of the Reflected Signal

The reflected phasor is attenuated by several factors:

### Fresnel Reflection Coefficients

**File:** `lib/snr/fwd/snr_fwd_reflected.m` → `snr_fwd_fresnelcoeff`

- Same-sense (RHCP→RHCP) and cross-sense (RHCP→LHCP) coefficients
- Depend on grazing angle and surface permittivity (dielectric constant)
- At low elevation (grazing incidence), coefficients can be large (strong reflection)

### Surface Roughness

**File:** `lib/snr/fwd/snr_fwd_roughness.m` → `get_roughness`

- Reduces coherent reflection; effect increases with `height_std/λ` and steeper incidence
- At low elevation (grazing), roughness has less impact than at high elevation

### Antenna Gain Pattern

**File:** `lib/snr/fwd/snr_fwd_antenna.m`

- Direct path: typically near boresight (zenith) → high gain
- Reflected path: arrives from low elevation (near horizon) → often lower gain
- Antenna pattern thus attenuates the reflected signal relative to direct

### Reflected Phasor Assembly

**File:** `lib/snr/fwd/snr_fwd_reflected.m`

```
phasor_nongeom = antenna_gain * (incident * Fresnel_same + incident_cross * Fresnel_cross)
phasor_geom = phasor_delay * phasor_roughness * phasor_divergence
phasor_reflected = phasor_geom * phasor_nongeom
```

---

## 6. SNR Conversion

**File:** `lib/snr/fwd/snr_fwd_signal2snr.m`

```
CN0 = power_signal / density_noise
CPN = CN0 / bandwidth_noise
CPN_db = 10*log10(CPN) + tracking_losses + estimator_effects
```

---

## 7. Data Flow Summary

```
Input: sat.elev, sat.azim, height_ant, wavelength, surface properties, antenna pattern
    │
    ▼
snr_fwd_geometry
    → delay = 2*h*sin(elev)
    │
    ▼
snr_fwd_reflected
    → phasor_reflected = (Fresnel × roughness × antenna) × exp(j·phase_delay)
    │
    ▼
snr_fwd_direct
    → phasor_direct
    │
    ▼
snr_fwd_combine
    → phasor_composite = phasor_direct + phasor_reflected
    │
    ▼
power_composite = |phasor_composite|²
    │
    ▼
SNR_db = f(power_composite, noise, bandwidth, ...)
    │
    ▼
get_multipath_modulation(phasor_direct, phasor_reflected)
    → extracts oscillation amplitude
    │
    ▼
mplsqfourier(SNR_detrended, sin(elev), period=λ/(2h), ...)
    → fits sinusoid, retrieves height, amplitude, phase
```

---

## 8. Key Equations for Implementation

### Reflection delay (horizontal surface)
```
delay_m = 2 * h_m * sin(elev_deg * π/180)
```

### Interferometric phase
```
phase_rad = 2π * delay_m / wavelength_m
```

### Oscillation in sin(elevation)
```
freq = 2 * h / λ    # cycles per unit sin(elev)
SNR_oscillation ≈ A * cos(2π * freq * sin(elev) + φ)
```

### Fringe visibility (oscillation strength)
```
r = |reflected|² / |direct|²
visibility = 2√r / (1 + r)
```

---

## 9. Applying to Wireless InSite / Other Simulations

If you have forward scattering in Wireless InSite but **do not see multipath oscillations at low elevation**, possible causes and remedies:

### 1. **Missing coherent reflection**

Wireless InSite may be computing diffuse/scattered power rather than a coherent specular reflection. The oscillation requires **coherent** addition of direct + reflected phasors (amplitude and phase). Ensure:
- A specular reflection path is modeled (or added post-hoc)
- Phase is preserved (not just power)

### 2. **Add the interference model post-hoc**

If Wireless InSite gives you:
- Direct path power (or field)
- Reflected path power (or field)

You can add the interference yourself:
```
E_composite = E_direct + E_reflected * exp(j * 2π * delay/λ)
power_composite = |E_composite|²
```

The reflected field should include Fresnel coefficient, roughness, and antenna pattern if not already in the simulation.

### 3. **Use this model to correct your data**

- Run this forward model with your geometry (height, elevation, wavelength)
- Extract the expected oscillation via `get_multipath_modulation`
- Compare to your Wireless InSite or measured SNR
- Fit `mplsqfourier` to measured data to retrieve effective reflector height and validate

### 4. **Critical parameters**

- **Antenna height** `h` – drives oscillation frequency
- **Wavelength** `λ` – L1 ≈ 0.19 m, L2 ≈ 0.24 m
- **Elevation** – use sin(elev) as the independent variable for the oscillation
- **Fresnel coefficients** – depend on surface permittivity and grazing angle

---

## 10. File Reference

| Purpose | Path |
|---------|------|
| Forward simulation entry | `lib/snr/fwd/snr_fwd.m` |
| Reflection delay geometry | `lib/snr/util/get_reflection_delay.m` |
| Delay → phase phasor | `lib/snr/fwd/snr_fwd_delay.m` |
| Reflected signal (Fresnel, roughness, antenna) | `lib/snr/fwd/snr_fwd_reflected.m` |
| Combine direct + reflected | `lib/snr/fwd/snr_fwd_combine.m` |
| Oscillation extraction | `lib/snr/util/get_multipath_modulation.m` |
| Oscillation extraction (aux) | `lib/snr/util/get_multipath_modulation_aux.m` |
| Fringe visibility | `lib/snr/util/get_fringe_vis.m` |
| Fringe deviation (dB) | `lib/snr/util/get_fringe_dev.m` |
| Sinusoid fit (sin(elev)) | `lib/snr/util/mplsqfourier/mplsqfourier.m` |
| Fringe fitting | `lib/snr/bias/snr_bias_fit_fringe.m` |
| Sinusoid fit driver | `lib/snr/bias/snr_bias_fit_sinusoid.m` |
| SNR conversion | `lib/snr/fwd/snr_fwd_signal2snr.m` |
| Fresnel coefficients | `lib/external/fresnel.m` (and snr_fwd_fresnelcoeff) |
| Surface permittivity | `lib/snr/util/get_permittivity.m` |

---

## 11. References

- Nievinski, F.G. and Larson, K.M., "An open source GPS multipath simulator in Matlab/Octave," GPS Solutions (2014).
- Zavorotny, V.U. et al., "A Physical Model for GPS Multipath Caused by Land Reflections: Toward Bare Soil Moisture Retrievals," IEEE J-STARS (2010).
