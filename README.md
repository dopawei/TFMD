# TFMD: Time-Frequency Mode Decomposition

**Paper:** Zhou W., Li W.-J., Zhu D., Xu H., Ren W.-X. (2026). [Time–frequency mode decomposition for wind turbine vibration monitoring under variable speed operation](https://doi.org/10.1016/j.oceaneng.2026.128438). *Ocean Engineering*, 368, 128438.

**Preprint:** [arXiv:2507.11919](https://doi.org/10.48550/arXiv.2507.11919)

MATLAB implementation of **Time-Frequency Mode Decomposition (TFMD)** for multicomponent nonstationary signals, such as rotor-order components of wind turbine vibration under variable speed operation.


## What Is TFMD?

TFMD defines each mode as the time domain signal reconstructed by inverse STFT from one connected support region in the short-time Fourier transform (STFT) plane. Signal decomposition thereby becomes the estimation of an unknown number of connected support regions, so the number of modes is an output of the segmentation rather than an input.

**Results reported in the paper:**
- The number of modes is inferred from the segmentation and matched the ground truth in all six synthetic cases from 10 to 40 dB input SNR
- Taking the median over seven input SNR levels (10-40 dB), TFMD attained the lowest average mode error in all six synthetic cases, compared with EMD, VMD, ACMD, SET, and VGNMD
- The average runtime stayed below 0.06 s per case and was second shortest, after EMD, in Cases 2-6
- Applying TFMD again to the residual (residual decomposition) recovered the weaker 2P harmonics of all nine operating states in a laboratory wind turbine blade strain experiment

## Quick Start

```matlab
% Reproduce Tables 1 and 2 of the paper (six synthetic signals, noise free)
test;
```

Expected summary (relative errors in units of 10^-2):

| Case | N / N_f | E_rel,total | E_rel,avg |
|------|---------|-------------|-----------|
| 1 | 2 / 2 | 2.09 | 2.59 |
| 2 | 2 / 2 | 3.01 | 3.64 |
| 3 | 4 / 4 | 3.38 | 4.68 |
| 4 | 2 / 2 | 4.44 | 5.01 |
| 5 | 7 / 7 | 3.62 | 5.64 |
| 6 | 6 / 6 | 0.12 | 0.18 |

## Basic Usage

### Decompose a signal in 3 lines:

```matlab
fs = 1000;                          % Sampling frequency
signal = your_signal;               % Your signal data
[modes, reconstructed] = tfmd(signal, fs);
```

### With custom parameters:

```matlab
opts.G = 128;        % Gaussian window length (samples)
opts.alpha = 2.5;    % Gaussian window shape parameter
opts.rho = 0.90;     % Overlap ratio
opts.beta = 0.5;     % Dilation factor
opts.sigma = 1e-3;   % Minimum area ratio

[modes, reconstructed] = tfmd(signal, fs, opts);
```

### Residual decomposition

When weak components are masked by dominant ones, the same procedure can be applied to the residual of the first decomposition:

```matlab
[modes1, recon1] = tfmd(signal, fs, opts);            % first decomposition
[modes2, recon2] = tfmd(signal - recon1, fs, opts);   % residual decomposition
final = recon1 + recon2;
```

## Files

| File | Description |
|------|-------------|
| `tfmd.m` | Core TFMD algorithm |
| `generate_signal.m` | Six synthetic signals of the paper (Section 3.1) |
| `test.m` | Reproduces Tables 1 and 2 of the paper |

## Synthetic Signals

| Case | Signal | Duration | f_s | N |
|------|--------|----------|-----|---|
| 1 | Frequency separated chirps | 1 s | 1000 Hz | 2 |
| 2 | Sinusoidal FM components | 1 s | 1000 Hz | 2 |
| 3 | Four component mixture | 1 s | 1000 Hz | 4 |
| 4 | Low frequency chirp and AM tone | 1 s | 1000 Hz | 2 |
| 5 | Generalized nonlinear signal | 3 s | 1000 Hz | 7 |
| 6 | Synthetic signal with stepped operating states | 120 s | 50 Hz | 6 |

Case 6 contains six operating states of 20 s each, with piecewise constant fundamental frequencies of 1.20, 1.85, 2.55, 3.25, 3.95, and 4.65 Hz. `generate_signal(case_idx)` returns each case at the sampling frequency used in the paper.

## Method Overview

TFMD works in 6 steps:

1. **STFT** - Transform the signal to the time-frequency plane with a Gaussian window
2. **Coefficient selection** - Apply two-cluster k-means to the STFT magnitudes to select coefficients dominated by signal energy
3. **Connected component labeling** - Group adjacent selected coefficients into connected regions (8-connectivity)
4. **Filtering by region size** - Remove regions smaller than the minimum area ratio
5. **Mask dilation with conflict resolution** - Expand each retained region while leaving contested bins unassigned
6. **Inverse STFT** - Reconstruct one mode from each final mask

## Parameters

| Symbol | Option | Default | Description |
|--------|--------|---------|-------------|
| $G$ | `G` | 128 | Gaussian window length (samples) |
| $\alpha$ | `alpha` | 2.5 | Gaussian window shape parameter |
| $\rho$ | `rho` | 0.90 | Overlap ratio |
| $\beta$ | `beta` | 0.5 | Dilation factor |
| $\sigma$ | `sigma` | 1e-3 | Minimum area ratio |

The FFT size is max(256, 2^nextpow2(G)). The option names `window_length`, `overlap_ratio`, `expansion_factor`, and `min_pixel_ratio` are also accepted.

### When to Adjust

- **Closely spaced or low frequency components**: reduce `alpha` for a less strongly tapered window and finer frequency resolution. The laboratory blade strain signal in the paper (6 Hz sampling) used `alpha = 1.25`.
- **Low SNR**: smaller `beta` values gave lower mode errors at low input SNR.
- **Noise free signals**: larger `beta` values reduced the mode error in Cases 1-5. The default `beta = 0.5` is a compromise across noise conditions.
- **Weak components left in the residual**: apply residual decomposition (see above).

## Requirements

- MATLAB R2020a or later (tested on R2024b)
- Signal Processing Toolbox (for `stft`, `istft`, `chirp`, `tukeywin`)
- Image Processing Toolbox (for `bwlabel`, `imdilate`)
- Statistics and Machine Learning Toolbox (for `kmeans`)
- Wavelet Toolbox (for `wextend`)
- Communications Toolbox (for `awgn`, only when `test.m` is run with a finite `snr_db`)

## Data Availability

The synthetic signals are generated by `generate_signal.m`. The laboratory wind turbine blade strain data analyzed in the paper are not publicly available due to project confidentiality constraints.

## Example: Two-tone signal
```matlab
fs = 1000;
t = (0:1/fs:1-1/fs)';
signal = sin(2*pi*100*t) + 0.8*sin(2*pi*200*t);

[modes, recon] = tfmd(signal, fs);
fprintf('Found %d modes\n', length(modes));
% Output: Found 2 modes
```

## Citation

If you use this code, please cite:

```bibtex
@article{zhou2026tfmd,
  title   = {Time--frequency mode decomposition for wind turbine vibration monitoring under variable speed operation},
  author  = {Zhou, Wei and Li, Wei-Jian and Zhu, Desen and Xu, Hongbin and Ren, Wei-Xin},
  journal = {Ocean Engineering},
  volume  = {368},
  pages   = {128438},
  year    = {2026},
  doi     = {10.1016/j.oceaneng.2026.128438}
}
```

## Related Work
Zhou et al. (2022). Empirical Fourier decomposition. *Mech. Syst. Signal Process.*, 163, 108155.
