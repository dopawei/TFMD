function signal_data = generate_signal(case_idx, fs)
% GENERATE_SIGNAL Generate synthetic test signals
%
% Usage:
%   signal_data = generate_signal(case_idx, fs)
%
% Inputs:
%   case_idx - Signal case (1-6)
%   fs       - Sampling frequency (Hz)
%
% Outputs:
%   signal_data - Struct with fields:
%                 .clean: composite signal
%                 .components_gt: ground truth components (cell array)
%                 .t: time vector
%                 .fs: sampling frequency
%                 .name: signal name
%                 .num_gt: number of components

% Initialize
T_dur = 1.0;  % Default duration
t = [];
components = {};
sig = [];

switch case_idx
    case 1  % Frequency-Separated Chirps
        N = round(T_dur * fs);
        t = (0:N-1)' / fs;
        c1 = 1.0 * chirp(t, 20, t(end), 70, 'linear');
        c2 = 0.9 * chirp(t, 130, t(end), 180, 'quadratic');
        components = {c1, c2};
        sig = c1 + c2;
        name = 'Frequency-Separated Chirps';
        
    case 2  % Sinusoidal FM
        N = round(T_dur * fs);
        t = (0:N-1)' / fs;
        c1 = 1.2 * cos(2*pi*100*t + 15*sin(2*pi*2*t));
        c2 = 1.0 * cos(2*pi*250*t + 5*sin(2*pi*5*t));
        components = {c1, c2};
        sig = c1 + c2;
        name = 'Sinusoidal FM';
        
    case 3  % Four Components Mix
        N = round(T_dur * fs);
        t = (0:N-1)' / fs;
        c1 = 1.0 * chirp(t, 10, t(end), 40);
        c2 = 0.9 * sin(2*pi*100*t);
        idx3 = (t >= 0) & (t <= 0.7);
        c3 = zeros(N, 1);
        c3(idx3) = 1.1 * cos(2*pi*350*t(idx3) + 5*sin(2*pi*6*t(idx3)));
        idx4 = (t >= 0.6) & (t <= 0.9);
        c4 = zeros(N, 1);
        c4(idx4) = 1.2 * sin(2*pi*200*t(idx4)) .* tukeywin(sum(idx4), 0.25);
        components = {c1, c2, c3, c4};
        sig = c1 + c2 + c3 + c4;
        name = 'Four Components Mix';
        
    case 4  % Chirp and AM Tone
        N = round(T_dur * fs);
        t = (0:N-1)' / fs;
        c1 = 1.0 * chirp(t, 20, t(end), 80);
        c2 = 1.1 * (0.8 + 0.4*cos(2*pi*2*t)) .* sin(2*pi*200*t);
        components = {c1, c2};
        sig = c1 + c2;
        name = 'Chirp and AM Tone';
        
    case 5  % Generalized Nonlinear (7 components, 3 seconds)
        T_dur = 3.0;
        t1 = 0:1/fs:1.5; t1 = t1(1:end-1);
        t2 = 1.5:1/fs:3; t2 = t2(1:end-1);
        t = [t1, t2]';
        N = length(t);
        
        % Components 1-3
        c1 = cos(2*pi*(170*t + 20*t.^2 + 3*cos(3*pi*t)));
        c2 = zeros(N, 1);
        c2(1:length(t1)) = cos(2*pi*(75*t1 + 20*t1.^2));
        t3 = 1:1/fs:3; t3 = t3(1:end-1);
        c3 = zeros(N, 1);
        idx3 = find(t >= 1, 1);
        c3(idx3:end) = cos(2*pi*(10*t3 + 20*t3.^2 + 3*cos(3*pi*t3)));
        
        % Components 4-7 (dispersive)
        Nf = floor(N/2) + 1;
        c4 = ifft_dispersive(Nf, N, T_dur, 1/2, @(f) 30*exp(-1j*2*pi*(0.4*f + 2*cos(2*pi*f/100))));
        c5 = ifft_dispersive(Nf, N, T_dur, 3/5, @(f) 30*exp(-1j*2*pi*(0.8*f + 0.0005*f.^2)));
        c6 = ifft_dispersive(Nf, N, T_dur, 7/10, @(f) 30*exp(-1j*2*pi*(1.8*f + 2*cos(2*pi*f/100))));
        c7 = ifft_dispersive(Nf, N, T_dur, 8/10, @(f) 30*exp(-1j*2*pi*(2.2*f + 0.0005*f.^2)));
        
        components = {c1, c2, c3, c4, c5, c6, c7};
        sig = c1 + c2 + c3 + c4 + c5 + c6 + c7;
        name = 'Generalized Nonlinear';
        
    case 6  % Two Simple Tones
        T_dur = 10.0;
        N = round(T_dur * fs);
        t = (0:N-1)' / fs;
        c1 = 1.0 * sin(2*pi*100*t);
        c2 = 0.8 * sin(2*pi*200*t);
        components = {c1, c2};
        sig = c1 + c2;
        name = 'Two Simple Tones';
        
    otherwise
        error('Invalid case_idx: %d (must be 1-6)', case_idx);
end

% Package output
signal_data.clean = sig(:);
signal_data.components_gt = components;
signal_data.t = t(:);
signal_data.fs = fs;
signal_data.name = name;
signal_data.num_gt = length(components);
signal_data.T_dur = T_dur;

end

function comp = ifft_dispersive(Nf, N, T_dur, ratio, func)
% Helper for Case 5 dispersive components
    idx_start = floor(ratio * Nf);
    idx_end = Nf - 1;
    f = (idx_start:idx_end) / T_dur;
    spec_pos = [zeros(1, idx_start), func(f)];
    
    % Hermitian symmetry
    if mod(N, 2) == 0
        spec_full = [spec_pos, conj(fliplr(spec_pos(2:end-1)))];
    else
        spec_full = [spec_pos, conj(fliplr(spec_pos(2:end)))];
    end
    
    comp = real(ifft(spec_full, N));
    comp = comp(:);
end