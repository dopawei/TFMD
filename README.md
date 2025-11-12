# TFMD: Time-Frequency Mode Decomposition
**Paper:** [Time-Frequency Mode Decomposition: A Morphological Segmentation Framework for Signal Analysis](https://doi.org/10.48550/arXiv.2507.11919) (arXiv 2025)

MATLAB implementation of **Time-Frequency Mode Decomposition (TFMD)** - a morphological segmentation framework for automatic multicomponent signal analysis.
...

## What is TFMD?

TFMD decomposes complex signals into individual components by treating signal decomposition as an image segmentation problem in the time-frequency domain. It automatically identifies and separates signal modes without knowing the number of components beforehand.

**Key advantages:**
- Automatic component detection (no need to specify number of modes)
- Noise robust (SNR: 10-40 dB)
- Fast (2nd fastest among benchmark methods)
- High accuracy in mode reconstruction

## Quick Start

```matlab
% Run test suite on 6 synthetic signals
test;
```

## Basic Usage

### Decompose a signal in 3 lines:

```matlab
fs = 1000;                          % Sampling frequency
signal = your_signal;               % Your signal data
[modes, reconstructed] = tfmd(signal, fs);
```

### With custom parameters:

```matlab
opts.G = 128;        % Window length
opts.alpha = 2.5;    % Gaussian shape
opts.rho = 0.90;     % Overlap ratio
opts.beta = 0.5;     % Expansion factor
opts.sigma = 1e-3;   % Size filter threshold

[modes, reconstructed] = tfmd(signal, fs, opts);
```

## Files

| File | Description |
|------|-------------|
| `tfmd.m` | Core TFMD algorithm (~150 lines) |
| `generate_signal.m` | Generate 6 test signals (~120 lines) |
| `test.m` | Complete test suite (~130 lines) |

## Method Overview

TFMD works in 6 steps:

1. **STFT** - Transform to time-frequency domain
2. **K-means** - Identify high-energy regions
3. **CCL** - Label connected components
4. **Filter** - Remove noise artifacts
5. **ICD** - Expand masks competitively
6. **ISTFT** - Reconstruct individual modes

## Parameters

### Core Parameters (Manuscript Notation)

| Symbol | Parameter | Default | Description |
|--------|-----------|---------|-------------|
| $G$ | Window length | 128 | STFT window size (samples) |
| $\alpha$ | Shape parameter | 2.5 | Gaussian window shape |
| $\rho$ | Overlap ratio | 0.90 | STFT window overlap |
| $\beta$ | Expansion factor | 0.5 | ICD mask expansion |
| $\sigma$ | Size threshold | 1e-3 | Minimum component size |

### When to Adjust

- **Closely spaced frequencies**: Reduce `alpha` to 1.5-2.0
- **High noise**: Reduce `beta` to 0.2-0.3
- **Low noise**: Increase `beta` to 0.8-1.0

## Requirements

- MATLAB R2020a or later
- Signal Processing Toolbox (for `stft`, `istft`)
- Image Processing Toolbox (for `bwlabel`)
- Statistics Toolbox (for `kmeans`)

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

```bibtex
@article{zhou2025tfmd,
  title={Time-Frequency Mode Decomposition: A Morphological Segmentation Framework for Signal Analysis},
  author={Zhou, Wei and Li, Wei-Jian and Zhu, Desen and Xu, Hongbin and Ren, Wei-Xin},
  year={2025},
  eprint={2507.11919},
  archivePrefix={arXiv},
  primaryClass={eess.SP}
}
```

## Related Work:
Zhou et al. (2022). Empirical Fourier decomposition. *Mech. Syst. Signal Process.*, 163, 108155.
