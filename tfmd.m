function [components, reconstructed_signal] = tfmd(input_signal, fs, options)
% TFMD - Time-Frequency Mode Decomposition
%
% Usage:
%   [components, reconstructed_signal] = tfmd(signal, fs)
%   [components, reconstructed_signal] = tfmd(signal, fs, options)
%
% Inputs:
%   signal  - Input signal (vector)
%   fs      - Sampling frequency (Hz)
%   options - Optional parameters (struct):
%             .window_length: STFT window length (default: 128)
%             .alpha: Gaussian window shape (default: 2.5)
%             .overlap_ratio: Window overlap (default: 0.90)
%             .beta: Expansion factor (default: 0.5)
%             .sigma: Size filter threshold (default: 1e-3)
%
% Outputs:
%   components - Cell array of extracted modes
%   reconstructed_signal - Sum of all modes

% Default parameters (manuscript notation)
if nargin < 3 || isempty(options)
    options = struct();
end

% Accept both old and new parameter names for backward compatibility
G = get_param(options, {'G', 'window_length'}, 128);
alpha = get_param(options, {'alpha'}, 2.5);
rho = get_param(options, {'rho', 'overlap_ratio'}, 0.90);
beta = get_param(options, {'beta', 'expansion_factor'}, 0.5);
sigma = get_param(options, {'sigma', 'min_pixel_ratio'}, 1e-3);

% Helper function to get parameter with multiple possible names
function val = get_param(opts, names, default)
    val = default;
    for i = 1:length(names)
        if isfield(opts, names{i})
            val = opts.(names{i});
            return;
        end
    end
end

% Prepare signal
input_signal = input_signal(:);
N = length(input_signal);
pad_len = floor(G / 2);
padded = wextend(1, 'sym', input_signal, pad_len);

% STFT
H = round(G * (1 - rho));
K = max(256, 2^nextpow2(G));
win = gausswin(G, alpha);
[S, F, T] = stft(padded, fs, 'Window', win, 'OverlapLength', G-H, ...
                 'FFTLength', K, 'FrequencyRange', 'centered');

% Extract positive frequencies
[Kf, M] = size(S);
k_dc = floor(Kf/2) + 1;
S_pos = S(k_dc:end, :);
magnitude = abs(S_pos);

% K-means clustering to find signal cores
features = magnitude(:);
try
    labels = kmeans(features, 2, 'MaxIter', 100, 'Replicates', 1);
catch
    labels = kmeans(features, 2, 'MaxIter', 100);
end
labels = reshape(labels, size(magnitude));
means = grpstats(features, labels(:), 'mean');
[~, signal_id] = max(means);
mask_core = (labels == signal_id);

% Connected-component labeling
[L0, N0] = bwlabel(mask_core, 8);

% Size filtering
min_size = round(numel(L0) * sigma);
L = L0;
for i = 1:N0
    if sum(L(:) == i) < min_size
        L(L == i) = 0;
    end
end
[L, Nf] = bwlabel(L > 0, 8);

% Iterative Competitive Dilation (ICD)
if Nf > 0
    % Calculate expansion radius for each component
    r = zeros(1, Nf);
    for i = 1:Nf
        area = sum(L(:) == i);
        r(i) = max(1, round(beta * sqrt(area / pi)));
    end
    
    % Expand regions
    L_expand = L;
    blocked = false(size(L));
    se = strel('square', 3);
    
    for iter = 1:max(r)
        claims = zeros(size(L));
        for i = 1:Nf
            if iter <= r(i)
                mask_i = (L_expand == i);
                dilated = imdilate(mask_i, se);
                frontier = dilated & ~mask_i & (L_expand == 0) & ~blocked;
                
                if any(frontier(:))
                    contested = frontier & (claims > 0);
                    new_claims = frontier & (claims == 0);
                    claims(new_claims) = i;
                    claims(contested) = -1;
                end
            end
        end
        
        blocked(claims == -1) = true;
        accept = (claims > 0) & ~blocked;
        if ~any(accept(:))
            break;
        end
        L_expand(accept) = claims(accept);
    end
else
    L_expand = L;
end

% Reconstruct modes
components = cell(1, Nf);
for i = 1:Nf
    % Create mask with Hermitian symmetry
    mask = false(Kf, M);
    mask(k_dc:end, :) = (L_expand == i);
    if mod(K, 2) == 0
        mask(2:k_dc-1, :) = flipud(mask(k_dc+1:end, :));
    else
        mask(1:k_dc-1, :) = flipud(mask(k_dc+1:end, :));
    end
    
    % ISTFT
    S_masked = S .* mask;
    recon_pad = real(istft(S_masked, fs, 'Window', win, 'OverlapLength', G-H, ...
                           'FFTLength', K, 'FrequencyRange', 'centered'));
    
    % Remove padding
    if length(recon_pad) >= pad_len + N
        components{i} = recon_pad(pad_len+1 : pad_len+N);
    else
        components{i} = zeros(N, 1);
    end
end

% Total reconstruction
if Nf > 0
    reconstructed_signal = sum(cat(2, components{:}), 2);
else
    reconstructed_signal = zeros(N, 1);
end

fprintf('TFMD: Extracted %d modes from signal (N=%d, fs=%d Hz)\n', Nf, N, fs);

end