clear, clc, close all;

%% Removing signals from muscle movement:

% Load the ECG signal data
load('ecg.mat');

% Extract the ECG signal and adjust for amplification factor
EKG1 =  ecg ./ 500;

% Calculate time vector based on sampling rate
fs = 500; % Sampling frequency (Hz)
t = (0:length(ecg)-1) / fs; % Time vector

% Plot the original EKG signal
figure;
subplot(221)
plot(t, EKG1);
xlabel('Time (seconds)');
ylabel('Voltage (V)');
grid on;
title('Original EKG Signal measured @ leads');
xlim([0, 1]); % Zoom into one period of the signal

% Perform Fourier Transform
N = length(ecg);
f = linspace(-fs/2, fs/2, N); % Frequency axis for FFT

% Compute the FFT of the signal
ecg_fft = fftshift(fft(ecg));

subplot(222)
plot(f, real(ecg_fft));
xlabel('Frequency (Hz)');
ylabel('Voltage (V)');
grid on;
title('ECG Signal in frequency domain');


% Set frequencies below 0.5 Hz to zero (low-pass filtering)
ecg_fft(f < 0.5) = 0;

subplot(224)
plot(f, real(ecg_fft));
xlabel('Frequency (Hz)');
ylabel('Voltage (V)');
grid on;
title('Removing frequency < 0.5 Hz');

% Inverse FFT to obtain the filtered signal
ecg1 = real(ifft(ifftshift(ecg_fft)));

% Plot the filtered ECG signal
subplot(223)
plot(t, real(ecg1)); % Use real part to avoid complex artifacts
xlabel('Time (seconds)');
ylabel('Voltage (V)');
grid on;
title('ECG Signal with Muscle Signals Removed');
xlim([0, 1]); % Zoom into one period of the signal

%% Removing 50 Hz interference:

% Design a Notch filter to remove 50 Hz interference
f0 = 50; % Frequency to notch out (Hz)
w0 = f0 / (fs/2); % Normalized frequency
Q = 30; % Quality factor for the Notch filter

% Design the Notch filter using the second-order section (SOS) structure
[b, a] = iirnotch(w0, w0/Q);
figure;
subplot(311)
phasez(b, a, 1024, fs);
title('Response of Notch Filter');
xlabel('Frequency (Hz)');
ylabel('Phase (Radians)');
grid on;

% Apply the Notch filter to the ECG signal
ecg2 = filtfilt(b, a, ecg1); % Zero-phase filtering

% Plot the original ECG signal
subplot(312)
plot(t, ecg);
xlabel('Time (seconds)');
ylabel('Voltage (V)');
grid on;
title('Original EKG Signal measured @ leads');
xlim([0, 1]); % Zoom into one period of the signal

% Plot the filtered ECG signal
subplot(313)
plot(t, ecg2);
grid on;
xlabel('Time (seconds)');
ylabel('Voltage (V)');
title('ECG Signal after 50 Hz Notch Filtering');
xlim([0, 1]); % Zoom into one period of the signal


%% Increasing the signal-to-noise ratio:

% Define different cutoff frequencies for the low-pass filter
cutoff_freqs = [20, 30, 40, 50]; % Hz

figure;
% Plot the original ECG signal
subplot(321)
plot(t, ecg);
xlabel('Time (seconds)');
ylabel('Voltage (V)');
grid on;
title('Original ECG Signal');
xlim([0, 1]); % Zoom into one period of the signal

% Apply low-pass filtering with different cutoff frequencies
for i = 1:length(cutoff_freqs)
    % Design a Butterworth low-pass filter
    cutoff = cutoff_freqs(i);
    [b, a] = butter(4, cutoff / (fs/2), 'low'); % 4th-order Butterworth filter
    
    % Apply the low-pass filter to the ECG signal
    filtered_signal = filtfilt(b, a, ecg2);
    
    % Plot the filtered ECG signal with specific color
	subplot(3,2,i+2)
    plot(t, filtered_signal);
	xlabel('Time (seconds)');
	ylabel('Voltage (V)');
	grid on;
	title(['Cutoff = ', num2str(cutoff), 'Hz']);
	xlim([0, 1]); % Zoom into one period of the signal
	hold on;
end
hold off; % Disable hold on

cutoff = 30;
[b, a] = butter(4, cutoff / (fs/2), 'low'); % 4th-order Butterworth filter

% Apply the low-pass filter to the ECG signal
ecg3 = filtfilt(b, a, ecg2);

% Plot the filtered ECG signal with specific color
subplot(322)
plot(t, ecg3);
xlabel('Time (seconds)');
ylabel('Voltage (V)');
grid on;
title('A compromise on the cut-off frequency, we found = 30 Hz');
xlim([0, 1]); % Zoom into one period of the signal

%% Finding the heart rate using autocorrelation:

% Compute the cross-correlation and lags for the original ECG signal
[ECG_autoc_original, lags] = xcorr(ecg);

% Plot the autocorrelation of the original ECG signal
figure
subplot(2,1,1)
plot(lags/fs, ECG_autoc_original);
grid on;
title('Autocorrelation of Original ECG Signal');
xlabel('Lag');
ylabel('Autocorrelation');

% Compute the cross-correlation and lags for ECG3 signal
[ECG_autoc3, lags] = xcorr(ecg3);

% Plot the autocorrelation of ECG3 signal
subplot(2,1,2)
plot(lags/fs, ECG_autoc3);
grid on;
title('Autocorrelation of ECG3 Signal');
xlabel('Lag');
ylabel('Autocorrelation');

% Find the global maximum of autocorrelation for the original ECG signal
[ECG_autoc_original_gmax, ECG_original_gmax_loc] = max(ECG_autoc_original);

% Find local peaks after the global maximum for the original EKG1 signal
[peaks, locations] = findpeaks(ECG_autoc_original((ECG_original_gmax_loc+1):end));
local_max_original_indices = find(ECG_autoc_original == max(peaks));
local_max_original_index = max(local_max_original_indices) - ECG_original_gmax_loc;

% Compute heart rate (in bpm) based on the local max index and sampling frequency (Fs) for the original ECG signal
heart_rate_original = (60 * fs) / local_max_original_index;

% Find the global maximum of autocorrelation for ECG3 signal
[ECG3_gmax, ECG3_gmax_loc] = max(ECG_autoc3);

% Find local peaks after the global maximum for ECG3 signal
[peaks3, locations3] = findpeaks(ECG_autoc3((ECG3_gmax_loc+1):end));
local_max3_indices = find(ECG_autoc3 == max(peaks3));
local_max3_index = max(local_max3_indices) - ECG3_gmax_loc;

% Compute heart rate (in bpm) based on the local max index and sampling frequency (Fs) for ECG3 signal
heart_rate3 = (60 * fs) / local_max3_index;

% Display the estimated heart rate
disp(['Estimated Heart Rate of original ECG: ' num2str(heart_rate_original) ' bpm ,While the Estimated Heart Rate of ECG3: ' num2str(heart_rate3) ' bpm']);

%% Finding the QRS complex:
figure;
 [qrs_amp_raw,qrs_i_raw,delay,mean_RR] = pan_tompkin(ecg3,fs,1);
