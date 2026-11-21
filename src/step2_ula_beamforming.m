%% 26 GHz mmWave LEO SatCom — Step 2: ULA Beamforming & Pointing-Error Analysis
% Course: ELEC-H415 (Wireless Communications / Propagation Modeling)
%
% Evaluates lambda/2-spaced Uniform Linear Arrays (ULAs) of N = [4, 8, 16, 200]
% elements for both isotropic and microstrip patch elements U(theta) = sin^2(theta):
%   1. Steered polar gain patterns at theta_m = [90, 45, 30, 0] deg
%   2. Received power and Shannon capacity under 0 deg and 2 deg pointing error
%   3. Bonus comparison of total pass data volume across pointing errors
%      [0, 0.2, 1, 2] deg (isotropic vs. patch arrays)
%
% Outputs: 5 figures + console summary table (1000-antenna batch sweep omitted).

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
tic;

%% Parameters
fc = 26e9;                 % Carrier frequency [Hz]
c = 3e8;                   % Speed of light [m/s]
lambda = c/fc;             % Wavelength [m]
beta = 2*pi/ lambda;       % Wavenumber [1/m]
omega = 2*pi*fc;           % Angular frequency [rad/s]
BW = 40e6;                 % Bandwidth [Hz]
EIRPW = 10e3;              % Transmitter EIRP [W] (eq 2.47)
EIRPdBW = 40;              % Transmitter EIRP [dBW]
EIRPdBm = 70;              % Transmitter EIRP [dBm]
F_db = 2;                  % Receiver Noise Figure [dB] (eq 3.92)
k = 1.379e-23;             % Boltzmann's constant [W Hz^-1 K^-1]
T_receiver = 275;          % Receiver noise temperature [K] (eq 3.92)
h_s = 500e3;               % Satellite orbit altitude [m]
T_rev = 5668;              % Satellite revolution period [s]
Z_0 = 377;                 % Free-space impedance [Ohm]
R_E = 6371e3;              % Earth radius [m]
theta = (0:0.5:180)*pi/180;% Elevation angle grid [rad]
d = sqrt( ((R_E+h_s)^2) - ( (R_E^2) * (cos(theta).^2) ) ) - (R_E.*sin(theta)); % Slant range [m]

% ULA parameters
antenna_space = lambda/2;
U_0 = 1;
size_array = [4; 8; 16; 200];              % Number of array elements
theta_m = [90, 45, 30, 0] * pi/180;        % Beam-steering angles [rad]

% Satellite pass timing weights
alpha = asin( (d/(R_E + h_s)) .* cos(theta) );
delta_alpha = alpha(1:end-1) - alpha(2:end);
delta_tn = T_rev * (delta_alpha/(2*pi));

%% Isotropic ULA: Steered Polar Gain Patterns
delta = (beta*antenna_space.*cos(theta_m))';
psi = beta*antenna_space.*cos(theta) - delta;
U_theta = U_0;

n = arrayfun(@(x) (0:(x-1)), size_array, 'UniformOutput', false);
n = cellfun(@(x) permute(x, [1 3 2]), n, 'UniformOutput', false);

Cell_Array_Factor = cellfun(@(x) sum(exp(1j * (x-1) .* psi), 3), n, 'UniformOutput', false);
Array_Factor = cat(3, Cell_Array_Factor{:});

U_array = arrayfun(@abs, Array_Factor).^2 .* U_theta;
P_rx = trapz(theta, U_array .* sin(theta), 2);
Gain_rx = U_array ./ ( (1/2) * P_rx );
Gain_reshaped_rx = permute(Gain_rx, [2 1 3]);

figure;
for idx = 1:numel(size_array)
    subplot(2, 2, idx);
    polarplot(theta, squeeze(Gain_reshaped_rx(:, :, idx)));
    title(sprintf('Polar plot %d isotropic antennas', size_array(idx)));
    thetalim([0 180]);
end
sgtitle('Isotropic ULA Radiation Patterns (\theta_m = 90^\circ, 45^\circ, 30^\circ, 0^\circ)');

%% Isotropic ULA: Received Power & Capacity with Pointing Error ([0, 0.2, 1, 2] deg)
pointing_error_deg = [0, 0.2, 1, 2];
pointing_error = pointing_error_deg * pi/180;
pointing_error = permute(pointing_error, [1 4 3 2]);
delta = beta*antenna_space.*cos(theta' + pointing_error);
psi = beta*antenna_space.*cos(theta) - delta;

Cell_Array_Factor = cellfun(@(x) sum(exp(1j * (x-1) .* psi), 3), n, 'UniformOutput', false);
Array_Factor = cat(3, Cell_Array_Factor{:});
U_array = arrayfun(@abs, Array_Factor).^2 .* U_theta;
sin_theta = sin(theta);
P_rx = trapz(theta, U_array .* sin_theta, 2);

[row, ~] = size(U_array);
for i = 1:row
    U_array(i, 1, :, :) = U_array(i, i, :, :);
end
U_array = U_array(:, 1, :, :);
Gain_rx = U_array ./ ( (1/2) * P_rx );
Gain_rx_iso = permute(Gain_rx, [3 1 4 2]);

P_rxW = EIRPW .* Gain_rx_iso .* (lambda / (4*pi) )^2 .* (1./d).^2;
P_rx_dBW = 10*log10( P_rxW );
P_rx_dBm = 10*log10( P_rxW / (1e-3) );

SNR_db = P_rx_dBW - F_db - 10*log10(k*T_receiver*BW);
SNR_lin = 10.^(SNR_db/10);
C = BW * log2(1+SNR_lin);
transmitted_data_iso = sum(((C(:,1:end-1,:) + C(:,2:end,:))/2) .* delta_tn, 2);
transmitted_data_iso = squeeze(transmitted_data_iso);

C_reshaped = permute(C, [2 1 3]);
P_rx_reshaped = permute(P_rx_dBm, [2 1 3]);
theta_deg = theta * 180 / pi;

figure;
subplot(2, 2, 1);
plot(theta_deg, P_rx_reshaped(:, :, 1)); grid on;
title('Received Power — Isotropic ULA (0^\circ Error)');
xlabel('Elevation angle \theta [deg]'); ylabel('P_{rx} [dBm]');
legend('4 antennas', '8 antennas', '16 antennas', '200 antennas', 'Location', 'south');

subplot(2, 2, 2);
plot(theta_deg, C_reshaped(:, :, 1) / 1e6); grid on;
title('Channel Capacity — Isotropic ULA (0^\circ Error)');
xlabel('Elevation angle \theta [deg]'); ylabel('Capacity [Mbit/s]');
legend('4 antennas', '8 antennas', '16 antennas', '200 antennas', 'Location', 'northwest');

subplot(2, 2, 3);
plot(theta_deg, P_rx_reshaped(:, :, 4)); grid on;
title('Received Power — Isotropic ULA (2^\circ Error)');
xlabel('Elevation angle \theta [deg]'); ylabel('P_{rx} [dBm]');
legend('4 antennas', '8 antennas', '16 antennas', '200 antennas', 'Location', 'south');

subplot(2, 2, 4);
plot(theta_deg, C_reshaped(:, :, 4) / 1e6); grid on;
title('Channel Capacity — Isotropic ULA (2^\circ Error)');
xlabel('Elevation angle \theta [deg]'); ylabel('Capacity [Mbit/s]');
legend('4 antennas', '8 antennas', '16 antennas', '200 antennas', 'Location', 'northwest');

%% Patch ULA: Steered Polar Gain Patterns (U_theta = U_0 sin^2(theta))
U_theta = U_0 * sin(theta).^2;
delta = (beta*antenna_space.*cos(theta_m))';
psi = beta*antenna_space.*cos(theta) - delta;

Cell_Array_Factor = cellfun(@(x) sum(exp(1j * x .* psi), 3), n, 'UniformOutput', false);
Array_Factor = cat(3, Cell_Array_Factor{:});
U_array = arrayfun(@abs, Array_Factor).^2 .* U_theta;
P_rx = trapz(theta, U_array .* sin(theta), 2);
Gain_rx = U_array ./ ( (1/2) * P_rx );
Gain_reshaped_rx = permute(Gain_rx, [2, 1, 3]);

figure;
for idx = 1:numel(size_array)
    subplot(2, 2, idx);
    polarplot(theta, squeeze(Gain_reshaped_rx(:, :, idx)));
    title(sprintf('Polar plot %d patch antennas', size_array(idx)));
    thetalim([0 180]);
end
sgtitle('Patch ULA Radiation Patterns (\theta_m = 90^\circ, 45^\circ, 30^\circ, 0^\circ)');

%% Patch ULA: Received Power & Capacity with Pointing Error ([0, 0.2, 1, 2] deg)
delta = beta*antenna_space.*cos((theta)' + pointing_error);
psi = beta*antenna_space.*cos(theta) - delta;

Cell_Array_Factor = cellfun(@(x) sum(exp(1j * (x-1) .* psi), 3), n, 'UniformOutput', false);
Array_Factor = cat(3, Cell_Array_Factor{:});
U_array = arrayfun(@abs, Array_Factor).^2 .* U_theta;
P_rx = trapz(theta, U_array .* sin(theta), 2);

[row, ~] = size(U_array);
for i = 1:row
    U_array(i, 1, :, :) = U_array(i, i, :, :);
end
U_array = U_array(:, 1, :, :);
Gain_rx = U_array ./ ( (1/2) * P_rx );
Gain_rx_patch = permute(Gain_rx, [3 1 4 2]);

P_rxW = EIRPW .* Gain_rx_patch .* (lambda / (4*pi) )^2 .* (1./d).^2;
P_rxW(:, length(theta), :) = 0;
P_rx_dBW = 10*log10( P_rxW );
P_rx_dBm = 10*log10( P_rxW / (1e-3) );

SNR_db = P_rx_dBW - F_db - 10*log10(k*T_receiver*BW);
SNR_lin = 10.^(SNR_db/10);
C = BW * log2(1+SNR_lin);
transmitted_data_patch = sum(((C(:,1:end-1,:) + C(:,2:end,:))/2) .* delta_tn, 2);
transmitted_data_patch = squeeze(transmitted_data_patch);

C_reshaped = permute(C, [2 1 3]);
P_rx_reshaped = permute(P_rx_dBm, [2 1 3]);

figure;
subplot(2, 2, 1);
plot(theta_deg, P_rx_reshaped(:, :, 1)); grid on;
title('Received Power — Patch ULA (0^\circ Error)');
xlabel('Elevation angle \theta [deg]'); ylabel('P_{rx} [dBm]');
legend('4 antennas', '8 antennas', '16 antennas', '200 antennas', 'Location', 'south');

subplot(2, 2, 2);
plot(theta_deg, C_reshaped(:, :, 1) / 1e6); grid on;
title('Channel Capacity — Patch ULA (0^\circ Error)');
xlabel('Elevation angle \theta [deg]'); ylabel('Capacity [Mbit/s]');
legend('4 antennas', '8 antennas', '16 antennas', '200 antennas', 'Location', 'northwest');

subplot(2, 2, 3);
plot(theta_deg, P_rx_reshaped(:, :, 4)); grid on;
title('Received Power — Patch ULA (2^\circ Error)');
xlabel('Elevation angle \theta [deg]'); ylabel('P_{rx} [dBm]');
legend('4 antennas', '8 antennas', '16 antennas', '200 antennas', 'Location', 'south');

subplot(2, 2, 4);
plot(theta_deg, C_reshaped(:, :, 4) / 1e6); grid on;
title('Channel Capacity — Patch ULA (2^\circ Error)');
xlabel('Elevation angle \theta [deg]'); ylabel('Capacity [Mbit/s]');
legend('4 antennas', '8 antennas', '16 antennas', '200 antennas', 'Location', 'northwest');

%% Bonus Comparison: Total Pass Data vs. Pointing Error
figure;
subplot(1, 2, 1);
plot(pointing_error_deg, transmitted_data_iso.' / 1e9, '-o', 'LineWidth', 1.5); grid on;
xlabel('Pointing error [deg]'); ylabel('Total transmitted data [Gbit]');
title('Isotropic ULA: Pointing Error Impact');
legend('4 antennas', '8 antennas', '16 antennas', '200 antennas', 'Location', 'northeast');

subplot(1, 2, 2);
plot(pointing_error_deg, transmitted_data_patch.' / 1e9, '-s', 'LineWidth', 1.5); grid on;
xlabel('Pointing error [deg]'); ylabel('Total transmitted data [Gbit]');
title('Patch ULA: Pointing Error Impact');
legend('4 antennas', '8 antennas', '16 antennas', '200 antennas', 'Location', 'northeast');
sgtitle('Bonus: Beamforming Sensitivity to Satellite Tracking Error');

fprintf('Step 2 Summary — Total Transmitted Data [Gbit] across [0, 0.2, 1, 2] deg:\n');
for idx = 1:numel(size_array)
    fprintf('  N=%3d | Iso: %s | Patch: %s\n', size_array(idx), ...
        sprintf('%6.2f ', transmitted_data_iso(idx,:)/1e9), ...
        sprintf('%6.2f ', transmitted_data_patch(idx,:)/1e9));
end
toc;
