%% Simplified TFMD Test Suite - All 6 Signal Cases
% Focused evaluation of TFMD on all synthetic signals

clear; close all; clc;

%% Global Parameters
fs = 1000;  % Sampling frequency (Hz)
num_cases = 6;  % Total number of test cases

fprintf('=== Simplified TFMD Test Suite ===\n');
fprintf('Testing all %d synthetic signal cases\n', num_cases);

%% TFMD Parameters (Manuscript Specifications)
base_options = struct();
base_options.window_length = 128;          % L_w = 128 samples
base_options.win_type = 'gaussian';        % Gaussian analysis window
base_options.alpha = 2.5;                  % Shape parameter α = 2.5
base_options.overlap_ratio = 115/128;      % 115 samples overlap → 89.8%
base_options.threshold_factor = 2.0;       % C_thresh = 2
base_options.min_component_size = 10;      % P_abs = 10 pixels
base_options.min_component_ratio = 0.005;  % P_rel = 0.005
base_options.denoise_filter_size = [3, 3]; % U × V = 3 × 3 kernel

fprintf('TFMD Parameters:\n');
fprintf('  - Gaussian window: L_w = %d samples, α = %.1f\n', base_options.window_length, base_options.alpha);
fprintf('  - Overlap: 115 samples (%.1f%%)\n', base_options.overlap_ratio*100);
fprintf('  - Threshold: C_thresh = %.1f\n', base_options.threshold_factor);
fprintf('  - Pixel thresholds: P_abs = %d, P_rel = %.3f\n', base_options.min_component_size, base_options.min_component_ratio);
fprintf('----------------------------------------\n\n');

%% Storage for Results
results = struct();

%% Main Test Loop
for case_idx = 1:num_cases
    fprintf('=== Testing Case %d ===\n', case_idx);
    
    %% Generate Signal
    signal_data = generate_signal(case_idx, fs);
    
    fprintf('Signal: %s\n', signal_data.name);
    fprintf('Ground truth components: %d\n', signal_data.num_gt);
    
    %% Set Case-Specific Parameters
    tfmd_options = base_options;
    
    % Special parameter for Case 5
    if case_idx == 5
        tfmd_options.alpha = 2.0;  % α = 2.0 for Case 5
        fprintf('Note: Using α = 2.0 for Case 5\n');
    end
    
    %% Apply TFMD
    fprintf('Applying TFMD...\n');
    [components, reconstructed_signal] = tfmd(signal_data.clean, fs, tfmd_options);
    
    %% Calculate Performance Metrics
    N_gt = signal_data.num_gt;
    N_f = length(components);
    
    % Total reconstruction error (E_rel,total)
    error_total = norm(signal_data.clean - reconstructed_signal) / norm(signal_data.clean);
    
    % Individual component matching and errors using minimum error method
    if N_f > 0 && N_gt > 0
        [component_errors, ~] = match_components_min_error(signal_data.components_gt, components);
        error_avg_mode = mean(component_errors);
    else
        component_errors = [];
        error_avg_mode = NaN;
    end
    
    %% Store Results
    results(case_idx).case_name = signal_data.name;
    results(case_idx).N_gt = N_gt;
    results(case_idx).N_f = N_f;
    results(case_idx).error_total = error_total;
    results(case_idx).error_avg_mode = error_avg_mode;
    results(case_idx).signal_data = signal_data;
    results(case_idx).components = components;
    results(case_idx).reconstructed_signal = reconstructed_signal;
    
    %% Display Case Results
    fprintf('Results:\n');
    fprintf('  Components: GT=%d, TFMD=%d\n', N_gt, N_f);
    fprintf('  Total error (E_rel,total): %.6f\n', error_total);
    if ~isnan(error_avg_mode)
        fprintf('  Avg mode error (E_rel,avg): %.6f\n', error_avg_mode);
    end
    fprintf('--------------------\n\n');
    
    %% Create Individual Case Figure
    create_case_figure(case_idx, results(case_idx));
end

%% Generate Summary Table
fprintf('=== TFMD Performance Summary ===\n');
fprintf('%-25s | %3s | %3s | %12s | %12s\n', ...
    'Case Name', 'N_g', 'N_f', 'E_rel,total', 'E_rel,avg');
fprintf('%s\n', repmat('-', 1, 70));

for i = 1:num_cases
    fprintf('%-25s | %3d | %3d | %12.2e | %12.2e\n', ...
        results(i).case_name, results(i).N_gt, results(i).N_f, ...
        results(i).error_total, results(i).error_avg_mode);
end
fprintf('%s\n', repmat('-', 1, 70));

%% Performance Statistics
total_components_gt = sum([results.N_gt]);
total_components_found = sum([results.N_f]);
detection_rate = total_components_found / total_components_gt * 100;
mean_total_error = mean([results.error_total]);

fprintf('\nOverall Performance Statistics:\n');
fprintf('  Total GT components: %d\n', total_components_gt);
fprintf('  Total TFMD components: %d\n', total_components_found);
fprintf('  Overall detection rate: %.1f%%\n', detection_rate);
fprintf('  Mean total error: %.6f\n', mean_total_error);

fprintf('\n=== TFMD Test Suite Completed ===\n');

%% Helper Functions

function [component_errors, matched_pairs] = match_components_min_error(gt_components, tfmd_components)
    % Simple component matching using minimum L2 error method
    
    N_gt = length(gt_components);
    N_tfmd = length(tfmd_components);
    
    if N_gt == 0 || N_tfmd == 0
        component_errors = [];
        matched_pairs = [];
        return;
    end
    
    % Compute error matrix (L2 normalized error between all pairs)
    error_matrix = zeros(N_gt, N_tfmd);
    for i = 1:N_gt
        for j = 1:N_tfmd
            % Ensure same length by padding or truncating
            len_min = min(length(gt_components{i}), length(tfmd_components{j}));
            gt_comp = gt_components{i}(1:len_min);
            tfmd_comp = tfmd_components{j}(1:len_min);
            
            % Compute normalized L2 error
            error_matrix(i,j) = norm(gt_comp - tfmd_comp) / (norm(gt_comp) + eps);
        end
    end
    
    % Greedy matching: find minimum error pairs
    matched_pairs = [];
    component_errors = [];
    used_tfmd = false(1, N_tfmd);
    used_gt = false(1, N_gt);
    
    for match_idx = 1:min(N_gt, N_tfmd)
        % Find minimum error among unused components
        temp_error = error_matrix;
        temp_error(used_gt, :) = inf;
        temp_error(:, used_tfmd) = inf;
        
        [min_error, linear_idx] = min(temp_error(:));
        [gt_idx, tfmd_idx] = ind2sub(size(temp_error), linear_idx);
        
        matched_pairs(end+1,:) = [gt_idx, tfmd_idx];
        component_errors(end+1) = min_error;
        
        used_gt(gt_idx) = true;
        used_tfmd(tfmd_idx) = true;
    end
end

function create_case_figure(case_idx, result)
    % Create visualization for individual test case
    
    signal_data = result.signal_data;
    components = result.components;
    t = signal_data.t;
    N_gt = result.N_gt;
    N_f = result.N_f;
    
    % Determine optimal subplot layout based on total components
    total_components = max(N_gt, N_f);
    if total_components <= 2
        subplot_rows = 2; subplot_cols = 3;
    elseif total_components <= 4
        subplot_rows = 3; subplot_cols = 3;
    else
        subplot_rows = 3; subplot_cols = 4;
    end
    
    figure('Name', sprintf('Case %d: %s', case_idx, result.case_name), ...
           'Position', [50 + case_idx*30, 50 + case_idx*30, 1200, 800]);
    
    % Get component matching
    if N_f > 0 && N_gt > 0
        [component_errors, matched_pairs] = match_components_min_error(signal_data.components_gt, components);
    else
        component_errors = [];
        matched_pairs = [];
    end
    
    % 1. Original vs Reconstructed Signal
    subplot(subplot_rows, subplot_cols, [1, 2]);
    plot(t*1000, signal_data.clean, 'k-', 'LineWidth', 2, 'DisplayName', 'Original');
    hold on;
    plot(t*1000, result.reconstructed_signal, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Reconstructed');
    title(sprintf('Case %d: %s', case_idx, result.case_name), 'FontWeight', 'bold');
    xlabel('Time (ms)'); ylabel('Amplitude');
    legend('Location', 'best'); grid on; axis tight;
    
    % Add performance metrics
    text(0.02, 0.98, sprintf('GT:%d, Found:%d, Error:%.4f', N_gt, N_f, result.error_total), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'BackgroundColor', 'white');
    
    % 2. Reconstruction Error
    subplot(subplot_rows, subplot_cols, 3);
    error_signal = signal_data.clean - result.reconstructed_signal;
    plot(t*1000, error_signal, 'g-', 'LineWidth', 1.5);
    title(sprintf('Error (RMS=%.4f)', rms(error_signal)));
    xlabel('Time (ms)'); ylabel('Error'); grid on; axis tight;
    
    % 3. Component Comparisons
    max_display = min(max(N_gt, N_f), subplot_rows * subplot_cols - 3);
    
    % Display matched pairs first
    comp_idx = 0;
    for i = 1:size(matched_pairs, 1)
        if comp_idx >= max_display; break; end
        comp_idx = comp_idx + 1;
        
        subplot(subplot_rows, subplot_cols, 3 + comp_idx);
        gt_idx = matched_pairs(i, 1);
        tfmd_idx = matched_pairs(i, 2);
        
        plot(t*1000, signal_data.components_gt{gt_idx}, 'b-', 'LineWidth', 2, 'DisplayName', 'GT');
        hold on;
        plot(t*1000, components{tfmd_idx}, 'r--', 'LineWidth', 1.5, 'DisplayName', 'TFMD');
        
        title(sprintf('GT%d ↔ TFMD%d (E=%.3f)', gt_idx, tfmd_idx, component_errors(i)));
        xlabel('Time (ms)'); ylabel('Amplitude');
        if comp_idx <= 2; legend('Location', 'best', 'FontSize', 8); end
        grid on; axis tight;
    end
    
    % Display unmatched GT components
    unmatched_gt = setdiff(1:N_gt, matched_pairs(:,1));
    for i = 1:length(unmatched_gt)
        if comp_idx >= max_display; break; end
        comp_idx = comp_idx + 1;
        
        subplot(subplot_rows, subplot_cols, 3 + comp_idx);
        plot(t*1000, signal_data.components_gt{unmatched_gt(i)}, 'b-', 'LineWidth', 2);
        title(sprintf('Unmatched GT%d', unmatched_gt(i)));
        xlabel('Time (ms)'); ylabel('Amplitude'); grid on; axis tight;
    end
    
    % Display unmatched TFMD components
    unmatched_tfmd = setdiff(1:N_f, matched_pairs(:,2));
    for i = 1:length(unmatched_tfmd)
        if comp_idx >= max_display; break; end
        comp_idx = comp_idx + 1;
        
        subplot(subplot_rows, subplot_cols, 3 + comp_idx);
        plot(t*1000, components{unmatched_tfmd(i)}, 'r--', 'LineWidth', 2);
        title(sprintf('Unmatched TFMD%d', unmatched_tfmd(i)));
        xlabel('Time (ms)'); ylabel('Amplitude'); grid on; axis tight;
    end
    
    sgtitle(sprintf('TFMD Analysis - Case %d: %s', case_idx, result.case_name), 'FontWeight', 'bold');
end