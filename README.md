# TFMD: Time-Frequency Mode Decomposition

This repository contains the MATLAB implementation for the **Time-Frequency Mode Decomposition (TFMD)** method.

## Introduction

TFMD is a novel signal processing framework that reframes signal decomposition as a morphological segmentation task in the time-frequency domain. The method is designed to automatically decompose multicomponent, non-stationary signals into their constituent modes by identifying morphologically distinct regions in the time-frequency plane. TFMD leverages unsupervised clustering, connected-component labeling, and iterative competitive dilation to construct precise time-frequency masks for each component, enabling accurate mode extraction without requiring prior specification of the number of components.

## Methodology Overview

The TFMD pipeline consists of the following key stages:

1. **Time-Frequency Representation**: Transform the input signal into the time-frequency domain using the Short-Time Fourier Transform (STFT) with a Gaussian window.

2. **Core Region Identification**: Apply unsupervised k-means clustering (C=2) to identify high-energy core regions of signal components, separating them from background noise.

3. **Connected-Component Labeling**: Segment the identified core regions using connected-component labeling (CCL) with 8-connectivity to isolate individual component cores.

4. **Size Filtering**: Apply a size filter to remove spurious detections, retaining only regions that exceed a minimum area threshold.

5. **Iterative Competitive Dilation (ICD)**: Expand each core region using the ICD algorithm to construct complete time-frequency masks while preventing overlap between different modes.

6. **Mode Reconstruction**: Apply the constructed masks to the STFT coefficients and reconstruct individual modes via inverse STFT (ISTFT) with weighted overlap-add (WOLA) synthesis.

7. **Optional Two-Stage Refinement (TFMD²)**: For signals with components of varying energy levels, apply TFMD to the residual signal to extract weaker modes that may be close to the noise floor.

## Key Features

- **Automatic Component Detection**: Determines the number of components without prior specification
- **Noise Robustness**: Maintains stable performance across SNR levels from 10 dB to 40 dB
- **Computational Efficiency**: Second fastest among benchmark methods, orders of magnitude faster than advanced alternatives
- **High Decomposition Fidelity**: Superior individual mode reconstruction accuracy
- **Adaptive Segmentation**: Morphological approach naturally handles time-localized and non-stationary components
- **Two-Stage Capability**: TFMD² approach effectively extracts both dominant and weaker harmonic components

## File Structure

- `tfmd.m`: The core function implementing the Time-Frequency Mode Decomposition algorithm
- `generate_signal.m`: Script to generate the six synthetic test signals used in validation
- `test.m`: Comprehensive test suite demonstrating TFMD on various synthetic signals
  
## Requirements

- **MATLAB** R2020a or later
- **Signal Processing Toolbox™**: Required for `stft`, `istft`, and windowing functions
- **Image Processing Toolbox™**: Required for connected-component labeling (`bwlabel` or `bwconncomp`)
- **Statistics and Machine Learning Toolbox™**: Required for k-means clustering

## How to Use

### Quick Start: Run the Test Suite

The easiest way to see TFMD in action is to run the complete test suite on six diverse synthetic signals:

```matlab
% Run comprehensive validation on all test cases
test;
```

This will generate signals, apply TFMD, calculate performance metrics, and create visualization figures for each case.

### Basic Usage

For custom applications, follow these steps:

#### Step 1: Prepare Your Signal

```matlab
fs = 1000;  % Sampling frequency in Hz
% Load or generate your signal
x = your_signal_data;
```

#### Step 2: Set TFMD Parameters

```matlab
% Define TFMD parameters (optional - defaults provided)
tfmd_params = struct();
tfmd_params.G = 128;           % Window length (samples)
tfmd_params.alpha = 2.5;       % Gaussian window shape parameter
tfmd_params.rho = 0.9;         % Overlap ratio
tfmd_params.sigma = 1e-3;      % Minimum area ratio (for size filter)
tfmd_params.beta = 0.5;        % Expansion factor (for ICD)
```

**Parameter Guidelines**:
- `G`: Window length controls time-frequency resolution trade-off (default: 128)
- `alpha`: Shape parameter for Gaussian window; reduce to 1.5-2.0 for closely spaced frequencies (default: 2.5)
- `rho`: Overlap ratio for STFT; higher values improve synthesis quality (default: 0.9)
- `sigma`: Minimum area threshold removes spurious noise detections (default: 1e-3)
- `beta`: Controls expansion radius in ICD; 0.5 balances energy capture and noise robustness (default: 0.5)

#### Step 3: Apply TFMD

```matlab
% Single-stage decomposition
[modes, x_reconstructed, N_f] = tfmd(x, fs, tfmd_params);

% Display results
fprintf('TFMD identified %d modes.\n', N_f);
```

#### Step 4: Two-Stage Refinement (Optional)

For signals with components having vastly different energy levels:

```matlab
% First stage
[modes_stage1, x_recon1, N_f1] = tfmd(x, fs, tfmd_params);

% Compute residual
residual = x - x_recon1;

% Second stage on residual
[modes_stage2, x_recon2, N_f2] = tfmd(residual, fs, tfmd_params);

% Combine results
all_modes = [modes_stage1, modes_stage2];
x_final = x_recon1 + x_recon2;

fprintf('Total modes extracted: %d (Stage 1: %d, Stage 2: %d)\n', ...
        N_f1 + N_f2, N_f1, N_f2);
```

## Performance Characteristics

Based on comprehensive numerical validation:

- **Baseline Accuracy**: Total reconstruction error εrel,total < 0.045 across all test cases
- **Noise Robustness**: Maintains εrel,avg < 0.1 at SNR ≥ 20 dB
- **Automatic Detection**: 100% accuracy in determining the correct number of components
- **Computational Speed**: 
  - 13× faster than ACMD for complex signals
  - 152× faster than SET
  - 577× faster than VGNMD

## Citation

If you use this code in your research, please cite:

```bibtex
@article{zhou2025tfmd,
  title={Time-Frequency Mode Decomposition: A Morphological Segmentation Framework for Signal Analysis and Its Application},
  author={Zhou, Wei and Li, Wei-Jian and Zhu, Desen and Xu, Hongbin and Ren, Wei-Xin},
      year={2025},
      eprint={2507.11919},
      archivePrefix={arXiv},
      primaryClass={eess.SP}
}
```

## Related Publications

For theoretical background and detailed methodology, please refer to:
- Zhou, W., et al. (2022). Empirical Fourier decomposition: An accurate signal decomposition method for nonlinear and non-stationary time series analysis. *Mechanical Systems and Signal Processing*, 163, 108155.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or collaboration inquiries:
- **Wei Zhou**: zhouw6@szu.edu.cn

College of Civil and Transportation Engineering  
Shenzhen University, Shenzhen 518061, China

