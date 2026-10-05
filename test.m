%% TFMD Test Suite - Six synthetic signal cases (Experiment 1 of the paper)
% Reproduces Tables 1 and 2 of:
%   Zhou W, Li W-J, Zhu D, Xu H, Ren W-X. Time-frequency mode decomposition
%   for wind turbine vibration monitoring under variable speed operation.
%   Ocean Engineering 368 (2026) 128438. doi:10.1016/j.oceaneng.2026.128438
clear; close all; clc;

fprintf('====================================\n');
fprintf('        TFMD Test Suite\n');
fprintf('====================================\n\n');

%% Parameters
snr_db = Inf;  % Input SNR in dB (Inf = noise free, as in Tables 1 and 2)

% TFMD parameters (Section 3.3 of the paper)
options = struct();
options.G = 128;            % Gaussian window length (samples)
options.alpha = 2.5;        % Gaussian window shape parameter
options.rho = 0.90;         % Overlap ratio
options.beta = 0.5;         % Dilation factor
options.sigma = 1e-3;       % Minimum area ratio

fprintf('TFMD Parameters:\n');
fprintf('  G     = %d samples (window length)\n', options.G);
fprintf('  alpha = %.1f (Gaussian shape)\n', options.alpha);
fprintf('  rho   = %.2f (overlap ratio)\n', options.rho);
fprintf('  beta  = %.2f (dilation factor)\n', options.beta);
fprintf('  sigma = %.0e (minimum area ratio)\n\n', options.sigma);

% Published values (Tables 1 and 2, noise free, x 1e-2)
paper_err_total = [2.09 3.01 3.38 4.44 3.62 0.12];
paper_err_avg   = [2.59 3.64 4.68 5.01 5.64 0.18];

%% Test all cases
results = cell(1, 6);

for case_idx = 1:6
    fprintf('--- Case %d ---\n', case_idx);

    % Generate signal (each case uses its own sampling frequency)
    data = generate_signal(case_idx);
    fs = data.fs;
    fprintf('Signal: %s (N=%d, fs=%g Hz)\n', data.name, data.num_gt, fs);

    % Add noise
    rng(2026);
    if isfinite(snr_db)
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
    [err_modes, pairing] = match_and_evaluate(data.components_gt, modes);
    err_avg = mean(err_modes);

    % Store results
    results{case_idx} = struct('name', data.name, 'N_gt', data.num_gt, ...
        'N_f', N_f, 'err_total', err_total, 'err_avg', err_avg, ...
        'err_modes', err_modes, 'time', time);

    fprintf('Found: %d modes | E_total: %.2e | E_avg: %.2e | Time: %.3fs\n\n', ...
            N_f, err_total, err_avg, time);

    % Plot
    plot_results(case_idx, data, modes, recon, pairing, err_modes);
end

%% Summary table
fprintf('=====================================================================\n');
fprintf('                              Summary\n');
fprintf('=====================================================================\n');
fprintf('Case | N/N_f | E_total (x1e-2) | E_avg (x1e-2) | Paper E_total / E_avg | Time\n');
fprintf('-----|-------|-----------------|---------------|-----------------------|------\n');
for i = 1:6
    r = results{i};
    fprintf('%4d | %d/%d   | %15.2f | %13.2f | %10.2f / %-10.2f| %.3f\n', ...
            i, r.N_gt, r.N_f, 100*r.err_total, 100*r.err_avg, ...
            paper_err_total(i), paper_err_avg(i), r.time);
end
fprintf('=====================================================================\n\n');

fprintf('Per-mode relative errors E_rel,i (x1e-2):\n');
for i = 1:6
    fprintf('  Case %d: %s\n', i, sprintf('%.2f ', 100*results{i}.err_modes));
end
fprintf('\n');

%% Helper functions

function [errors, pairing] = match_and_evaluate(gt, modes)
    % Pairing rule of the paper (Section 3.2): each ground truth component,
    % in order, is assigned the still-unused mode with the largest absolute
    % Pearson correlation. Unpaired components count as error 1.
    N_gt = length(gt);
    N_modes = length(modes);
    errors = ones(1, N_gt);
    pairing = zeros(1, N_gt);

    C = zeros(N_gt, N_modes);
    for i = 1:N_gt
        for j = 1:N_modes
            len = min(length(gt{i}), length(modes{j}));
            if std(gt{i}(1:len)) > eps && std(modes{j}(1:len)) > eps
                r = corrcoef(gt{i}(1:len), modes{j}(1:len));
                C(i, j) = abs(r(1, 2));
            end
        end
    end

    used = false(1, N_modes);
    for i = 1:N_gt
        [~, order] = sort(C(i, :), 'descend');
        j = order(find(~used(order), 1, 'first'));
        if isempty(j)
            continue;
        end
        used(j) = true;
        pairing(i) = j;
        len = min(length(gt{i}), length(modes{j}));
        errors(i) = norm(gt{i}(1:len) - modes{j}(1:len)) / norm(gt{i}(1:len));
    end
end

function plot_results(case_idx, data, modes, recon, pairing, err_modes)
    % Original vs reconstruction, error, and each paired mode
    figure('Position', [50+case_idx*30, 50+case_idx*30, 1200, 800]);

    t = data.t;  % seconds
    N_gt = data.num_gt;
    N_f = length(modes);

    % Determine layout
    if max(N_gt, N_f) <= 3
        rows = 2; cols = 3;
    elseif max(N_gt, N_f) <= 6
        rows = 3; cols = 3;
    else
        rows = 3; cols = 4;
    end

    % 1. Original vs reconstructed
    subplot(rows, cols, [1 2]);
    plot(t, data.clean, 'k-', 'LineWidth', 2); hold on;
    plot(t, recon, 'r--', 'LineWidth', 1.5);
    title(sprintf('Case %d: %s', case_idx, data.name));
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('Original', 'Reconstructed');
    grid on;

    % 2. Error
    subplot(rows, cols, 3);
    plot(t, data.clean - recon, 'g-');
    title(sprintf('Error (RMS=%.4f)', rms(data.clean - recon)));
    xlabel('Time (s)'); grid on;

    % 3. Paired components
    plot_idx = 4;
    for i = 1:N_gt
        if plot_idx > rows*cols; break; end
        subplot(rows, cols, plot_idx);
        plot(t, data.components_gt{i}, 'b-', 'LineWidth', 2); hold on;
        if pairing(i) > 0
            plot(t, modes{pairing(i)}, 'r--', 'LineWidth', 1.5);
            title(sprintf('GT%d <-> Mode%d (E=%.2f%%)', i, pairing(i), 100*err_modes(i)));
            legend('GT', 'TFMD', 'Location', 'best');
        else
            title(sprintf('Unmatched GT%d', i), 'Color', 'r');
        end
        xlabel('Time (s)'); ylabel('Amplitude');
        grid on;
        plot_idx = plot_idx + 1;
    end

    % 4. Unpaired TFMD modes
    unpaired = setdiff(1:N_f, pairing(pairing > 0));
    for j = unpaired
        if plot_idx > rows*cols; break; end
        subplot(rows, cols, plot_idx);
        plot(t, modes{j}, 'r--', 'LineWidth', 2);
        title(sprintf('Unmatched Mode%d', j), 'Color', 'r');
        xlabel('Time (s)'); grid on;
        plot_idx = plot_idx + 1;
    end

    sgtitle(sprintf('Case %d: %s', case_idx, data.name));
end
