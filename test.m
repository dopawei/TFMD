%% TFMD Test Suite - All 6 Signal Cases
clear; close all; clc;

fprintf('====================================\n');
fprintf('        TFMD Test Suite\n');
fprintf('====================================\n\n');

%% Parameters
fs = 1000;  % Sampling frequency
snr_db = 100;  % SNR (set to 100 for near noise-free)

% TFMD parameters (manuscript notation)
options = struct();
options.G = 128;            % Window length (Eq. 2)
options.alpha = 2.5;        % Gaussian shape parameter α (Eq. 2)
options.rho = 115/128;      % Overlap ratio ρ (Eq. 3)
options.beta = 0.5;         % Expansion factor β (Eq. 7)
options.sigma = 1e-3;       % Minimum area ratio σ (Eq. 9)

fprintf('TFMD Parameters (Manuscript Notation):\n');
fprintf('  G     = %d samples (window length)\n', options.G);
fprintf('  alpha = %.1f (Gaussian shape)\n', options.alpha);
fprintf('  rho   = %.4f (overlap ratio: 115/128)\n', options.rho);
fprintf('  beta  = %.2f (expansion factor)\n', options.beta);
fprintf('  sigma = %.0e (min area ratio)\n\n', options.sigma);

%% Test all cases
results = cell(1, 6);

for case_idx = 1:6
    fprintf('--- Case %d ---\n', case_idx);
    
    % Generate signal
    data = generate_signal(case_idx, fs);
    fprintf('Signal: %s (N_gt=%d)\n', data.name, data.num_gt);
    
    % Add noise
    rng(133);
    if snr_db < 100
        signal = awgn(data.clean, snr_db, 'measured');
    else
        signal = data.clean;
    end
    
    % Run TFMD
    tic;
    [modes, recon] = tfmd(signal, fs, options);
    time = toc;
    
    % Evaluate
    N_f = length(modes);
    err_total = norm(data.clean - recon) / norm(data.clean);
    
    % Match components
    if N_f > 0 && data.num_gt > 0
        err_modes = match_and_evaluate(data.components_gt, modes);
        err_avg = mean(err_modes);
    else
        err_avg = NaN;
    end
    
    % Store results
    results{case_idx} = struct('name', data.name, 'N_gt', data.num_gt, ...
        'N_f', N_f, 'err_total', err_total, 'err_avg', err_avg, ...
        'time', time, 'data', data, 'modes', {modes}, 'recon', recon);
    
    fprintf('Found: %d modes | Error: %.4f | Time: %.3fs\n\n', ...
            N_f, err_total, time);
    
    % Plot
    plot_results(case_idx, data, modes, recon);
end

%% Summary table
fprintf('====================================\n');
fprintf('           Summary\n');
fprintf('====================================\n');
fprintf('Case | Name                    | GT | Found | Error    | Time\n');
fprintf('-----|-------------------------|-------|-------|----------|------\n');
for i = 1:6
    r = results{i};
    fprintf('%4d | %-23s | %3d | %5d | %.2e | %.3f\n', ...
            i, r.name, r.N_gt, r.N_f, r.err_total, r.time);
end
fprintf('====================================\n\n');

% Overall stats
all_gt = sum(cellfun(@(x) x.N_gt, results));
all_found = sum(cellfun(@(x) x.N_f, results));
fprintf('Total: %d/%d components detected (%.1f%%)\n', ...
        all_found, all_gt, 100*all_found/all_gt);

%% Helper functions

function errors = match_and_evaluate(gt, modes)
    % Simple greedy matching by minimum error
    N_gt = length(gt);
    N_modes = length(modes);
    errors = zeros(1, min(N_gt, N_modes));
    
    % Error matrix
    E = zeros(N_gt, N_modes);
    for i = 1:N_gt
        for j = 1:N_modes
            len = min(length(gt{i}), length(modes{j}));
            E(i,j) = norm(gt{i}(1:len) - modes{j}(1:len)) / norm(gt{i}(1:len));
        end
    end
    
    % Greedy matching
    used_gt = false(1, N_gt);
    used_modes = false(1, N_modes);
    for k = 1:min(N_gt, N_modes)
        E_temp = E;
        E_temp(used_gt, :) = inf;
        E_temp(:, used_modes) = inf;
        [err, idx] = min(E_temp(:));
        [i, j] = ind2sub(size(E), idx);
        errors(k) = err;
        used_gt(i) = true;
        used_modes(j) = true;
    end
end

function plot_results(case_idx, data, modes, recon)
    % Simple visualization with proper component matching
    figure('Position', [50+case_idx*30, 50+case_idx*30, 1200, 800]);
    
    t = data.t * 1000;  % to ms
    N_gt = data.num_gt;
    N_f = length(modes);
    
    % Determine layout
    if max(N_gt, N_f) <= 2
        rows = 2; cols = 3;
    elseif max(N_gt, N_f) <= 4
        rows = 3; cols = 3;
    else
        rows = 3; cols = 4;
    end
    
    % 1. Original vs reconstructed
    subplot(rows, cols, [1 2]);
    plot(t, data.clean, 'k-', 'LineWidth', 2); hold on;
    plot(t, recon, 'r--', 'LineWidth', 1.5);
    title(sprintf('Case %d: %s', case_idx, data.name));
    xlabel('Time (ms)'); ylabel('Amplitude');
    legend('Original', 'Reconstructed');
    grid on;
    
    % 2. Error
    subplot(rows, cols, 3);
    plot(t, data.clean - recon, 'g-');
    title(sprintf('Error (RMS=%.4f)', rms(data.clean - recon)));
    xlabel('Time (ms)'); grid on;
    
    % 3. Match components first
    matched_pairs = [];
    if N_f > 0 && N_gt > 0
        % Build error matrix
        E = zeros(N_gt, N_f);
        for i = 1:N_gt
            for j = 1:N_f
                len = min(length(data.components_gt{i}), length(modes{j}));
                E(i,j) = norm(data.components_gt{i}(1:len) - modes{j}(1:len)) / ...
                         norm(data.components_gt{i}(1:len));
            end
        end
        
        % Greedy matching
        used_gt = false(1, N_gt);
        used_modes = false(1, N_f);
        matched_pairs = zeros(min(N_gt, N_f), 3);  % [gt_idx, mode_idx, error]
        
        for k = 1:min(N_gt, N_f)
            E_temp = E;
            E_temp(used_gt, :) = inf;
            E_temp(:, used_modes) = inf;
            [err, idx] = min(E_temp(:));
            [i, j] = ind2sub(size(E), idx);
            matched_pairs(k, :) = [i, j, err];
            used_gt(i) = true;
            used_modes(j) = true;
        end
    end
    
    % 4. Plot matched components
    plot_idx = 4;
    num_matched = size(matched_pairs, 1);
    
    for k = 1:min(num_matched, rows*cols-3)
        subplot(rows, cols, plot_idx);
        
        gt_idx = matched_pairs(k, 1);
        mode_idx = matched_pairs(k, 2);
        error = matched_pairs(k, 3);
        
        plot(t, data.components_gt{gt_idx}, 'b-', 'LineWidth', 2); hold on;
        plot(t, modes{mode_idx}, 'r--', 'LineWidth', 1.5);
        
        title(sprintf('GT%d <-> Mode%d (E=%.3f)', gt_idx, mode_idx, error*100));
        xlabel('Time (ms)'); 
        ylabel('Amplitude');
        legend('GT', 'TFMD', 'Location', 'best');
        grid on;
        
        plot_idx = plot_idx + 1;
    end
    
    % 5. Plot unmatched GT components
    if N_gt > 0 && num_matched > 0
        unmatched_gt = setdiff(1:N_gt, matched_pairs(:,1));
        for i = 1:length(unmatched_gt)
            if plot_idx > rows*cols; break; end
            
            subplot(rows, cols, plot_idx);
            plot(t, data.components_gt{unmatched_gt(i)}, 'b-', 'LineWidth', 2);
            title(sprintf('Unmatched GT%d', unmatched_gt(i)), 'Color', 'r');
            xlabel('Time (ms)'); grid on;
            plot_idx = plot_idx + 1;
        end
    end
    
    % 6. Plot unmatched TFMD modes
    if N_f > 0 && num_matched > 0
        unmatched_modes = setdiff(1:N_f, matched_pairs(:,2));
        for i = 1:length(unmatched_modes)
            if plot_idx > rows*cols; break; end
            
            subplot(rows, cols, plot_idx);
            plot(t, modes{unmatched_modes(i)}, 'r--', 'LineWidth', 2);
            title(sprintf('Unmatched Mode%d', unmatched_modes(i)), 'Color', 'r');
            xlabel('Time (ms)'); grid on;
            plot_idx = plot_idx + 1;
        end
    end
    
    sgtitle(sprintf('Case %d: %s', case_idx, data.name));
end