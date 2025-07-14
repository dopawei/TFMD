function signal_data = generate_signal(signal_choice, fs)
% GENERATE_SIGNAL Generates synthetic signals as defined in the TFMD paper
%   
% Signal Cases (matching manuscript equations):
%       1: Frequency-Separated Chirps (Equation 15)
%       2: Sinusoidal FM Components (Equation 16) 
%       3: Four Components Mix (Equation 17)
%       4: Low-Frequency Chirp and AM Tone (Equation 18)
%       5: Generalized Nonlinear Signal (Equations 19-23)
%       6: Two Simple Tones (Equation 24)
%
% Inputs:
%       signal_choice: Integer index of the signal to generate (1-6)
%       fs: Sampling frequency in Hz
%
% Outputs:
%       signal_data: Structure containing:
%           .clean: Clean composite signal (column vector)
%           .components_gt: Ground truth components (cell array)
%           .t: Time vector (column vector)
%           .fs: Sampling frequency
%           .name: Signal name
%           .num_gt: Number of ground truth components
%           .T_dur: Signal duration

fprintf('--- Signal Generation ---\n');

switch signal_choice
    case 1 % Case 1: Frequency-Separated Chirps (Equation 15)
        fprintf('Generating Signal Case 1: Frequency-Separated Chirps (fs=%.0f Hz)\n', fs);
        T_dur = 1.0;
        t = (0:1/fs:T_dur-1/fs)';
        
        % Linear chirp: 20 Hz to 70 Hz
        x_1_1 = 1.0 * cos(2*pi * (20*t + 25*t.^2));
        
        % Quadratic chirp: 130 Hz to 205 Hz  
        x_1_2 = 0.9 * cos(2*pi * (130*t + 25*t.^3));
        
        signal_clean = x_1_1 + x_1_2;
        components_gt = {x_1_1, x_1_2};
        signal_name = 'Frequency-Separated Chirps';
        
    case 2 % Case 2: Sinusoidal FM Components (Equation 16)
        fprintf('Generating Signal Case 2: Sinusoidal FM Components (fs=%.0f Hz)\n', fs);
        T_dur = 1.0;
        t = (0:1/fs:T_dur-1/fs)';
        
        % Component 1: fc=100Hz, fm=2Hz, frequency deviation=30Hz
        x_2_1 = 1.2 * cos(2*pi * 100*t + 15*sin(2*pi * 2*t));
        
        % Component 2: fc=250Hz, fm=5Hz, frequency deviation=25Hz  
        x_2_2 = 1.0 * cos(2*pi * 250*t + 5*sin(2*pi * 5*t));
        
        signal_clean = x_2_1 + x_2_2;
        components_gt = {x_2_1, x_2_2};
        signal_name = 'Sinusoidal FM Components';
        
    case 3 % Case 3: Four Components Mix (Equation 17)
        fprintf('Generating Signal Case 3: Four Components Mix (fs=%.0f Hz)\n', fs);
        T_dur = 1.0;
        t = (0:1/fs:T_dur-1/fs)';
        
        % Component 1: Linear chirp 10-40 Hz
        x_3_1 = 1.0 * cos(2*pi * (10*t + 15*t.^2));
        
        % Component 2: Pure tone at 100 Hz
        x_3_2 = 0.9 * sin(2*pi * 100*t);
        
        % Component 3: Time-limited FM signal (0 to 0.7s)
        x_3_3 = zeros(size(t));
        idx_3 = (t >= 0) & (t <= 0.7);
        t_active = t(idx_3) - 0.1;  % Shifted time for the component
        t_active(t_active < 0) = 0;  % Ensure non-negative
        x_3_3(idx_3) = 1.1 * cos(2*pi * 350*t_active + 5*sin(2*pi * 6*t_active));
        
        % Component 4: Transient AM burst (0.6 to 0.9s)
        x_3_4 = zeros(size(t));
        idx_4 = (t >= 0.6) & (t <= 0.9);
        N_burst = sum(idx_4);
        if N_burst > 0
            tukey_window = tukeywin(N_burst, 0.25);
            x_3_4(idx_4) = 1.2 * tukey_window .* sin(2*pi * 200*t(idx_4));
        end
        
        signal_clean = x_3_1 + x_3_2 + x_3_3 + x_3_4;
        components_gt = {x_3_1, x_3_2, x_3_3, x_3_4};
        signal_name = 'Four Components Mix';
        
    case 4 % Case 4: Low-Frequency Chirp and AM Tone (Equation 18)
        fprintf('Generating Signal Case 4: Low-Frequency Chirp and AM Tone (fs=%.0f Hz)\n', fs);
        T_dur = 1.0;
        t = (0:1/fs:T_dur-1/fs)';
        
        % Component 1: Linear chirp 20-80 Hz
        x_4_1 = 1.0 * cos(2*pi * (20*t + 30*t.^2));
        
        % Component 2: AM tone (200 Hz carrier, 2 Hz modulation)
        x_4_2 = 1.1 * (0.8 + 0.4*cos(2*pi * 2*t)) .* sin(2*pi * 200*t);
        
        signal_clean = x_4_1 + x_4_2;
        components_gt = {x_4_1, x_4_2};
        signal_name = 'Low-Frequency Chirp and AM Tone';
        
    case 5 % Case 5: Generalized Nonlinear Signal (Equations 19-23)
        fprintf('Generating Signal Case 5: Generalized Nonlinear Signal (fs=%.0f Hz)\n', fs);
        
        % Generate time vectors for two segments
        t1 = 0:1/fs:1.5; t1 = t1(1:end-1);
        t2 = 1.5:1/fs:3; t2 = t2(1:end-1);
        t = [t1 t2]; % Row vector
        t = t(:); % Convert to column vector
        Nt = length(t);
        T_dur = t(end) + 1/fs;
        
        % Additional time segments for component construction
        t11 = 0:1/fs:1; t11 = t11(1:end-1);
        t22 = 1:1/fs:3; t22 = t22(1:end-1);
        
        % Component 1: Full duration nonlinear chirp
        Sig1 = cos(2*pi*(170*t(:)'+20*t(:)'.^2+3*cos(3*pi*t(:)')));
        
        % Component 2: First half duration (0 to 1.5s)
        Sig21 = cos(2*pi*(75*t1+20*t1.^2)); 
        Sig22 = zeros(1,length(t2)); 
        Sig2 = [Sig21,Sig22];
        
        % Component 3: Second half duration (1 to 3s) using absolute time
        Sig31 = zeros(1,length(t11)); 
        Sig32 = cos(2*pi*(10*t22+20*t22.^2+3*cos(3*pi*t22))); 
        Sig3 = [Sig31,Sig32];
        
        % Frequency domain parameters for spectral components
        Nf = floor(Nt/2)+1; 
        T_fft = T_dur;
        
        % Component 4: Spectral definition (Equation 20)
        f41_idx = 0 : floor(1*Nf/2)-1; 
        f42_idx = floor(1*Nf/2) : Nf-1; 
        f42 = f42_idx / T_fft; 
        Ds41 = zeros(1,length(f41_idx)); 
        Ds42 = 30*exp(-1j*2*pi*(0.4*f42+2*cos(2*pi*f42/100))); 
        CDs4 = [complex(Ds41) Ds42];
        
        % Component 5: Spectral definition (Equation 21)
        f51_idx = 0 : floor(3*Nf/5)-1; 
        f52_idx = floor(3*Nf/5) : Nf-1; 
        f52 = f52_idx / T_fft; 
        Ds51 = zeros(1,length(f51_idx)); 
        Ds52 = 30*exp(-1j*2*pi*(0.8*f52+0.0005*f52.^2)); 
        CDs5 = [complex(Ds51) Ds52];
        
        % Component 6: Spectral definition (Equation 22)
        f61_idx = 0 : floor(7*Nf/10)-1; 
        f62_idx = floor(7*Nf/10) : Nf-1; 
        f62 = f62_idx / T_fft; 
        Ds61 = zeros(1,length(f61_idx)); 
        Ds62 = 30*exp(-1j*2*pi*(1.8*f62+2*cos(2*pi*f62/100))); 
        CDs6 = [complex(Ds61) Ds62];
        
        % Component 7: Spectral definition (Equation 23)
        f71_idx = 0 : floor(8*Nf/10)-1; 
        f72_idx = floor(8*Nf/10) : Nf-1; 
        f72 = f72_idx / T_fft; 
        Ds71 = zeros(1,length(f71_idx)); 
        Ds72 = 30*exp(-1j*2*pi*(2.2*f72+0.0005*f72.^2)); 
        CDs7 = [complex(Ds71) Ds72];
        
        % Time domain synthesis using inverse FFT with Hermitian symmetry
        ifftSig4 = real(ifft(hermitian_fft(CDs4, Nt), Nt, 2)); Sig4 = ifftSig4;
        ifftSig5 = real(ifft(hermitian_fft(CDs5, Nt), Nt, 2)); Sig5 = ifftSig5;
        ifftSig6 = real(ifft(hermitian_fft(CDs6, Nt), Nt, 2)); Sig6 = ifftSig6;
        ifftSig7 = real(ifft(hermitian_fft(CDs7, Nt), Nt, 2)); Sig7 = ifftSig7;
        
        % Combine all components
        Sig = Sig1+Sig2+Sig3+Sig4+Sig5+Sig6+Sig7;
        signal_clean = Sig(:);
        components_gt = {Sig1(:), Sig2(:), Sig3(:), Sig4(:), Sig5(:), Sig6(:), Sig7(:)};
        signal_name = 'Generalized Nonlinear Signal';
        
    case 6 % Case 6: Two Simple Tones (Equation 24)
        fprintf('Generating Signal Case 6: Two Simple Tones (fs=%.0f Hz)\n', fs);
        T_dur = 1.0;
        t = (0:1/fs:T_dur-1/fs)';
        
        % Component 1: 100 Hz tone
        x_6_1 = 1.0 * sin(2*pi * 100*t);
        
        % Component 2: 200 Hz tone
        x_6_2 = 0.8 * sin(2*pi * 200*t);
        
        signal_clean = x_6_1 + x_6_2;
        components_gt = {x_6_1, x_6_2};
        signal_name = 'Two Simple Tones';
        
    otherwise
        error('generate_signal:invalid_choice', 'Invalid signal_choice: %d. Select 1-6.', signal_choice);
end

% Ensure all outputs are column vectors
signal_clean = signal_clean(:);
t = t(:);
for i = 1:length(components_gt)
    components_gt{i} = components_gt{i}(:);
end

% Package output structure
signal_data = struct();
signal_data.clean = signal_clean;
signal_data.components_gt = components_gt;
signal_data.t = t;
signal_data.fs = fs;
signal_data.name = signal_name;
signal_data.num_gt = length(components_gt);
signal_data.T_dur = T_dur;

fprintf('Signal generation complete: "%s" (N=%d, Duration=%.3fs, Components=%d)\n', ...
    signal_data.name, length(signal_data.clean), signal_data.T_dur, signal_data.num_gt);
fprintf('---------------------------\n\n');

end

%% Helper function for Case 5 spectral synthesis
function full_fft = hermitian_fft(positive_freq_spectrum, N)
    % HERMITIAN_FFT Ensures Hermitian symmetry for real signal IFFT
    % Creates full FFT spectrum from positive frequencies only
    %
    % Inputs:
    %   positive_freq_spectrum: Complex spectrum for positive frequencies
    %   N: Total length of desired FFT
    %
    % Output:
    %   full_fft: Full spectrum with Hermitian symmetry
    
    len_pos = length(positive_freq_spectrum);
    
    if mod(N, 2) == 0 % Even N
        % Expects DC, Freqs 1 to N/2-1, Nyquist (N/2) -> len = N/2 + 1
        if len_pos ~= N/2 + 1
            warning('hermitian_fft: Even N length mismatch. Expected %d, got %d.', N/2+1, len_pos);
            % Basic correction attempt
            if len_pos > N/2+1
                positive_freq_spectrum = positive_freq_spectrum(1:N/2+1);
            else
                positive_freq_spectrum = [positive_freq_spectrum, zeros(1, N/2+1-len_pos)];
            end
        end
        % Create full spectrum with conjugate symmetry
        full_fft = [positive_freq_spectrum, conj(fliplr(positive_freq_spectrum(2:end-1)))];
        
    else % Odd N
        % Expects DC, Freqs 1 to (N-1)/2 -> len = (N-1)/2 + 1
        if len_pos ~= (N-1)/2 + 1
            warning('hermitian_fft: Odd N length mismatch. Expected %d, got %d.', (N-1)/2+1, len_pos);
            if len_pos > (N-1)/2+1
                positive_freq_spectrum = positive_freq_spectrum(1:(N-1)/2+1);
            else
                positive_freq_spectrum = [positive_freq_spectrum, zeros(1, (N-1)/2+1-len_pos)];
            end
        end
        % Create full spectrum with conjugate symmetry
        full_fft = [positive_freq_spectrum, conj(fliplr(positive_freq_spectrum(2:end)))];
    end
    
    % Final length verification and correction
    if length(full_fft) ~= N
        warning('hermitian_fft: Final length mismatch after construction. Expected %d, got %d.', N, length(full_fft));
        % Attempt simple padding/truncation
        if length(full_fft) > N
            full_fft = full_fft(1:N);
        else
            full_fft = [full_fft, zeros(1, N - length(full_fft))];
        end
    end
end