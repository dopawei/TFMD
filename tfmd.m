function [components, reconstructed_signal] = tfmd(input_signal, fs, options)
% tfmd - Time-Frequency Mode Decomposition
%
% Input parameters:
%   input_signal    - Input signal (column vector)
%   fs             - Sampling frequency (Hz)
%   options        - Optional parameter structure
%
% Output parameters:
%   components     - Decomposed modal components (cell array)
%   reconstructed_signal - Reconstructed signal (column vector)
%
% Features:
%   - Gaussian window only
%   - No padding support
%   - Optimized for synthetic signal analysis
%
% Author: Wei Zhou, College of Civil and Transportation Engineering, Shenzhen University, Shenzhen, 518061, China
% Date: July 2025

%% Default parameter settings
params = struct();
params.window_length = 128;           % Window length (samples)
params.alpha = 2.5;                   % Gaussian window shape parameter
params.overlap_ratio = 0.90;          % 90% overlap
params.threshold_factor = 2.0;        % Threshold factor (C_thresh)
params.min_component_size = 10;       % Minimum connected component size (P_abs)
params.min_component_ratio = 0.005;   % Minimum connected component ratio (P_rel)
params.denoise_filter_size = [3, 3];  % Spectrogram smoothing filter size

% User parameter override
if nargin >= 3 && isstruct(options)
    fields = fieldnames(options);
    for i = 1:length(fields)
        if isfield(params, fields{i})
            params.(fields{i}) = options.(fields{i});
        end
    end
end

% Ensure input is a column vector
if size(input_signal, 2) > 1
    input_signal = input_signal(:);
end

original_length = length(input_signal);

%% Step 1: Compute Short-Time Fourier Transform (STFT)
window_length = min(params.window_length, length(input_signal));
overlap_length = round(window_length * params.overlap_ratio);
nfft = max(256, 2^nextpow2(window_length));

% Create Gaussian window
if license('test', 'Signal_Toolbox')
    window = gausswin(window_length, params.alpha);
else
    % Manual Gaussian window creation
    t_win = -(window_length-1)/2:(window_length-1)/2;
    window = exp(-0.5 * (params.alpha * t_win / ((window_length-1)/2)).^2)';
end

% Compute STFT
[S, F, T] = stft(input_signal, fs, 'Window', window, ...
    'OverlapLength', overlap_length, 'FFTLength', nfft, ...
    'FrequencyRange', 'centered');

%% Step 2: Extract non-negative frequency components
[num_freq_bins, num_time_frames] = size(S);
dc_index = floor(num_freq_bins/2) + 1;
positive_freq_indices = dc_index:num_freq_bins;

S_positive = S(positive_freq_indices, :);
F_positive = F(positive_freq_indices);

% Compute magnitude spectrogram
magnitude_spectrum = abs(S_positive);

%% Step 3: Spectrogram smoothing (denoising)
magnitude_spectrum_smoothed = magnitude_spectrum;
if all(params.denoise_filter_size > 1)
    try
        % Apply 2D smoothing filter
        filter_kernel = ones(params.denoise_filter_size) / prod(params.denoise_filter_size);
        if exist('imfilter', 'file')
            magnitude_spectrum_smoothed = imfilter(magnitude_spectrum, filter_kernel, 'replicate');
        else
            % Manual convolution if Image Processing Toolbox not available
            magnitude_spectrum_smoothed = conv2(magnitude_spectrum, filter_kernel, 'same');
        end
    catch
        warning('Spectrogram smoothing failed. Using original magnitude spectrum.');
    end
end

%% Step 4: Adaptive threshold calculation
magnitude_values = magnitude_spectrum_smoothed(:);
max_magnitude = max(magnitude_values);
median_magnitude = median(magnitude_values);

% Geometric mean threshold (Equation 13 in paper)
if median_magnitude < eps
    if max_magnitude < eps
        threshold = 0;
    else
        mean_magnitude = mean(magnitude_values);
        threshold = sqrt(max_magnitude * max(mean_magnitude, eps) / params.threshold_factor);
    end
else
    threshold = sqrt(max_magnitude * median_magnitude / params.threshold_factor);
end

%% Step 5: Initial binary mask generation
binary_mask = magnitude_spectrum_smoothed > threshold;

%% Step 6: Connected-component labeling (CCL)
try
    [labeled_components, num_components] = bwlabel(binary_mask, 8);
catch
    % If Image Processing Toolbox is not available, use simple CCL algorithm
    [labeled_components, num_components] = simple_connected_components(binary_mask);
end

%% Step 7: Size-based filtering of candidate masks
min_pixels = max(params.min_component_size, ...
    floor(params.min_component_ratio * numel(binary_mask)));

valid_components = [];
for i = 1:num_components
    component_mask = (labeled_components == i);
    if sum(component_mask(:)) >= min_pixels
        valid_components(end+1) = i;
    end
end

%% Step 8: Construct full frequency domain masks and reconstruct signals
num_valid_components = length(valid_components);
components = cell(1, num_valid_components);

for i = 1:num_valid_components
    comp_id = valid_components(i);
    
    % Get positive frequency mask
    positive_mask = (labeled_components == comp_id);
    
    % Construct full frequency domain mask (including negative frequencies)
    full_mask = false(num_freq_bins, num_time_frames);
    full_mask(positive_freq_indices, :) = positive_mask;
    
    % Symmetric extension for negative frequencies (conjugate symmetry)
    for k = 1:(dc_index-1)
        symmetric_index = 2*dc_index - k;
        source_index = symmetric_index - dc_index + 1;
        if source_index >= 1 && source_index <= size(positive_mask, 1)
            full_mask(k, :) = positive_mask(source_index, :);
        end
    end
    
    % Apply mask to extract masked TF component
    masked_stft = S .* full_mask;
    
    % Inverse STFT (ISTFT) for time-domain mode reconstruction
    try
        reconstructed_component = istft(masked_stft, fs, 'Window', window, ...
            'OverlapLength', overlap_length, 'FFTLength', nfft, ...
            'FrequencyRange', 'centered', 'ConjugateSymmetric', true);
        
        % Adjust length to match original signal
        if length(reconstructed_component) > original_length
            reconstructed_component = reconstructed_component(1:original_length);
        elseif length(reconstructed_component) < original_length
            reconstructed_component(end+1:original_length) = 0;
        end
        
        components{i} = reconstructed_component(:);
    catch ME
        warning('Failed to reconstruct component %d: %s', i, ME.message);
        components{i} = zeros(original_length, 1);
    end
end

%% Step 9: Final signal synthesis
if num_valid_components > 0
    reconstructed_signal = sum(cat(2, components{:}), 2);
else
    reconstructed_signal = zeros(original_length, 1);
end

% Ensure output is column vector with correct length
reconstructed_signal = reconstructed_signal(:);
if length(reconstructed_signal) ~= original_length
    if length(reconstructed_signal) > original_length
        reconstructed_signal = reconstructed_signal(1:original_length);
    else
        temp_signal = zeros(original_length, 1);
        temp_signal(1:length(reconstructed_signal)) = reconstructed_signal;
        reconstructed_signal = temp_signal;
    end
end

end

%% Auxiliary function: Simple connected-component labeling algorithm
function [labeled, num_labels] = simple_connected_components(binary_image)
% Simple 8-connected component labeling algorithm (without Image Processing Toolbox)

[rows, cols] = size(binary_image);
labeled = zeros(rows, cols);
num_labels = 0;

% 8-connected neighborhood offsets
neighbors = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1];

for i = 1:rows
    for j = 1:cols
        if binary_image(i,j) && labeled(i,j) == 0
            % Found new connected component
            num_labels = num_labels + 1;
            
            % Use depth-first search to label connected component
            stack = [i, j];
            while ~isempty(stack)
                current_i = stack(end, 1);
                current_j = stack(end, 2);
                stack(end, :) = [];
                
                if current_i >= 1 && current_i <= rows && ...
                   current_j >= 1 && current_j <= cols && ...
                   binary_image(current_i, current_j) && labeled(current_i, current_j) == 0
                    
                    labeled(current_i, current_j) = num_labels;
                    
                    % Add neighbors to stack
                    for k = 1:size(neighbors, 1)
                        ni = current_i + neighbors(k, 1);
                        nj = current_j + neighbors(k, 2);
                        stack(end+1, :) = [ni, nj];
                    end
                end
            end
        end
    end
end

end