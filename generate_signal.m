function signal_data = generate_signal(case_idx, fs)
% GENERATE_SIGNAL Generate the six synthetic test signals of the paper
%
% Zhou W, Li W-J, Zhu D, Xu H, Ren W-X. Time-frequency mode decomposition
% for wind turbine vibration monitoring under variable speed operation.
% Ocean Engineering 368 (2026) 128438. doi:10.1016/j.oceaneng.2026.128438
%
% Usage:
%   signal_data = generate_signal(case_idx)
%
% Inputs:
%   case_idx - Signal case (1-6), see Section 3.1 of the paper
%   fs       - (optional, ignored) kept for backward compatibility. Each
%              case uses the sampling frequency stated in the paper:
%              1000 Hz for Cases 1-5 and 50 Hz for Case 6.
%
% Outputs:
%   signal_data - Struct with fields:
%                 .clean: composite signal
%                 .components_gt: ground truth components (cell array)
%                 .t: time vector
%                 .fs: sampling frequency
%                 .name: signal name
%                 .num_gt: number of components
%
% Cases:
%   1: Frequency separated chirps,          1 s @ 1000 Hz, N = 2
%   2: Sinusoidal FM components,            1 s @ 1000 Hz, N = 2
%   3: Four component mixture,              1 s @ 1000 Hz, N = 4
%   4: Low frequency chirp and AM tone,     1 s @ 1000 Hz, N = 2
%   5: Generalized nonlinear signal,        3 s @ 1000 Hz, N = 7
%   6: Synthetic signal with stepped
%      operating states,                  120 s @   50 Hz, N = 6

fs_paper = 1000;
if case_idx == 6
    fs_paper = 50;
end
if nargin >= 2 && ~isempty(fs) && fs ~= fs_paper
    warning('generate_signal:fsIgnored', ...
        'Case %d uses fs = %g Hz as in the paper; input fs = %g Hz is ignored.', ...
        case_idx, fs_paper, fs);
end
fs = fs_paper;

switch case_idx
    case 1  % Frequency separated chirps, Eq. (25)
        T_dur = 1.0;
        N = round(T_dur * fs);
        t = (0:N-1)' / fs;
        c1 = 1.0 * chirp(t, 20, t(end), 70);
        c2 = 0.9 * chirp(t, 130, t(end), 180, 'quadratic');
        components = {c1, c2};
        name = 'Frequency separated chirps';

    case 2  % Sinusoidal FM components, Eq. (26)
        T_dur = 1.0;
        N = round(T_dur * fs);
        t = (0:N-1)' / fs;
        c1 = 1.2 * cos(2*pi*100*t + (30/2) * sin(2*pi*2*t));
        c2 = 1.0 * cos(2*pi*250*t + (25/5) * sin(2*pi*5*t));
        components = {c1, c2};
        name = 'Sinusoidal FM components';

    case 3  % Four component mixture, Eq. (27)
        T_dur = 1.0;
        N = round(T_dur * fs);
        t = (0:N-1)' / fs;
        c1 = 1.0 * chirp(t, 10, t(end), 40);
        c2 = 0.9 * sin(2*pi*100*t);
        idx3 = (t >= 0) & (t <= 0.7);
        c3 = zeros(N, 1);
        c3(idx3) = 1.1 * cos(2*pi*350*t(idx3) + (30/6) * sin(2*pi*6*t(idx3)));
        idx4 = (t >= 0.6) & (t <= 0.9);
        c4 = zeros(N, 1);
        c4(idx4) = 1.2 * sin(2*pi*200*t(idx4)) .* tukeywin(sum(idx4), 0.25);
        components = {c1, c2, c3, c4};
        name = 'Four component mixture';

    case 4  % Low frequency chirp and AM tone, Eq. (28)
        T_dur = 1.0;
        N = round(T_dur * fs);
        t = (0:N-1)' / fs;
        c1 = 1.0 * chirp(t, 20, t(end), 80);
        c2 = 1.1 * (0.8 + 0.4*cos(2*pi*2*t)) .* sin(2*pi*200*t);
        components = {c1, c2};
        name = 'Low frequency chirp and AM tone';

    case 5  % Generalized nonlinear signal, Eqs. (29)-(30)
        t1 = 0:1/fs:1.5; t1 = t1(1:end-1);
        t2 = 1.5:1/fs:3; t2 = t2(1:end-1);
        t_row = [t1, t2];
        t = t_row(:);
        N = length(t);
        T_dur = t(end) + 1/fs;

        t11 = 0:1/fs:1; t11 = t11(1:end-1);
        t22 = 1:1/fs:3; t22 = t22(1:end-1);

        % Components 1-3
        c1 = cos(2*pi*(170*t_row + 20*t_row.^2 + 3*cos(3*pi*t_row)));
        c2 = [cos(2*pi*(75*t1 + 20*t1.^2)), zeros(1, length(t2))];
        c3 = [zeros(1, length(t11)), cos(2*pi*(10*t22 + 20*t22.^2 + 3*cos(3*pi*t22)))];

        % Components 4-7 from prescribed complex spectra G_i(f)
        Nf = floor(N/2) + 1;
        c4 = ifft_dispersive(Nf, N, T_dur, 1/2,  @(f) 30*exp(-1j*2*pi*(0.4*f + 2*cos(2*pi*f/100))));
        c5 = ifft_dispersive(Nf, N, T_dur, 3/5,  @(f) 30*exp(-1j*2*pi*(0.8*f + 0.0005*f.^2)));
        c6 = ifft_dispersive(Nf, N, T_dur, 7/10, @(f) 30*exp(-1j*2*pi*(1.8*f + 2*cos(2*pi*f/100))));
        c7 = ifft_dispersive(Nf, N, T_dur, 8/10, @(f) 30*exp(-1j*2*pi*(2.2*f + 0.0005*f.^2)));

        components = {c1, c2, c3, c4, c5, c6, c7};
        name = 'Generalized nonlinear signal';

    case 6  % Synthetic signal with stepped operating states, Eqs. (31)-(34)
        state_duration = 20;                                  % s per operating state
        f0_states = [1.20, 1.85, 2.55, 3.25, 3.95, 4.65];     % Hz
        T_dur = state_duration * numel(f0_states);            % 120 s
        N = round(T_dur * fs);
        t = (0:N-1)' / fs;

        state_idx = min(floor(t / state_duration) + 1, numel(f0_states));
        f0 = f0_states(state_idx);
        theta = 2*pi*cumsum(f0(:)) / fs;

        components = cell(1, numel(f0_states));
        for s = 1:numel(f0_states)
            idx = (state_idx == s);
            local_time = t(idx) - min(t(idx));
            envelope = (1.0 + 0.04*cos(2*pi*0.06*local_time + 0.35*s)) .* tukeywin(sum(idx), 0.30);
            c = zeros(N, 1);
            c(idx) = envelope .* cos(theta(idx));
            components{s} = c;
        end
        name = 'Synthetic signal with stepped operating states';

    otherwise
        error('Invalid case_idx: %d (must be 1-6)', case_idx);
end

% Package output
sig = zeros(length(t), 1);
for k = 1:numel(components)
    components{k} = components{k}(:);
    sig = sig + components{k};
end

signal_data.clean = sig;
signal_data.components_gt = components;
signal_data.t = t(:);
signal_data.fs = fs;
signal_data.name = name;
signal_data.num_gt = numel(components);
signal_data.T_dur = T_dur;

end

function comp = ifft_dispersive(Nf, N, T_dur, ratio, func)
% Helper for Case 5: real component from a one-sided complex spectrum that is
% zero below floor(ratio*Nf) and equal to func(f) above it.
    idx_start = floor(ratio * Nf);
    f = (idx_start:Nf-1) / T_dur;
    spec_pos = [complex(zeros(1, idx_start)), func(f)];

    % Hermitian symmetry
    if mod(N, 2) == 0
        spec_full = [spec_pos, conj(fliplr(spec_pos(2:end-1)))];
    else
        spec_full = [spec_pos, conj(fliplr(spec_pos(2:end)))];
    end

    comp = real(ifft(spec_full, N, 2));
    comp = comp(:);
end
