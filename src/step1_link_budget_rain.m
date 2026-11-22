%% 26 GHz mmWave LEO SatCom — Step 1: Free-Space Link Budget & Rain Attenuation
% Course: ELEC-H415 (Wireless Communications / Propagation Modeling)
%
% Models a 26 GHz LEO downlink (h_s = 500 km, T_rev = 5668 s, BW = 40 MHz)
% from horizon to horizon (elevation theta = 0..180 deg) with a single
% lossless isotropic receive antenna:
%   1. ITU-style rain attenuation across exceedance probabilities (0..0.001%)
%   2. Free-space path loss L0, Friis received power, SNR, and Shannon capacity
%   3. Integrated pass data volume and constant-bitrate optimization bonus
%
% Outputs: 6 figures + workspace metrics (transmitted_data, total_data).

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end

%% Parameters
% Assumptions:
% - Horizontally polarized wave matched to RX polarization (phi = 0, theta-only)
% - Atmospheric attenuation neglected except in the Rain Attenuation section

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
kb = 1.379e-23;            % Boltzmann's constant [W Hz^-1 K^-1]
T_receiver = 275;          % Receiver noise temperature [K] (eq 3.92)
h_s = 500e3;               % Satellite orbit altitude [m]
T_rev = 5668;              % Satellite revolution period [s]
Z_0 = 377;                 % Free-space impedance [Ohm]
R_E = 6371e3;              % Earth radius [m]
theta = (0:180)*pi/180;    % Satellite elevation angle [rad]
d = sqrt( ((R_E+h_s)^2) - ( (R_E^2) * (cos(theta).^2) ) ) - (R_E.*sin(theta));  % Slant range [m]

%% Bonus 1.1: Rain Attenuation
h_r = 2.400;                        % Rain height [km]
R_001 = 30;                         % Rainfall rate exceeded 0.01% [mm/h]
exceedance = [0; 1; 0.1; 0.01; 0.001];
d_s = h_r ./ (sin(theta));          % Slant-path length through rain
k = 0.15 + 0.003* (cos(theta).^2);
alpha_rain = (0.15 + 0.007* (cos(theta).^2) ) ./ k;
Gamma = k .* R_001.^alpha_rain;
L_r_001 = Gamma .* d_s;
L_r_001(1:21) = L_r_001(21);
L_r_001(161:end) = L_r_001(161);
Attenuation_list_L_r = L_r_001 .* (exceedance/0.01).^(-0.65 - ( 0.03 * log(exceedance) ) + (0.05* log(L_r_001) ) );

figure();
plot(theta*180/pi, Attenuation_list_L_r');
xlabel('Elevation angle \theta [deg]');
ylabel('Attenuation [dB]');
title('Rain Attenuation vs. Satellite Elevation');
legend('No Attenuation', 'Rain Attenuation: 1%', 'Rain Attenuation: 0.1%', 'Rain Attenuation: 0.01%', 'Rain Attenuation: 0.001%', 'Location', 'best');
grid on;

%% Free-Space Path Loss
G_rx = 1; % Lossless single isotropic receive antenna gain
L0 = 20*log10(d) + 20*log10(fc) - 147.5 + Attenuation_list_L_r; % Path loss [dB] (eq 3.12)

figure();
plot(theta*180/pi, L0');
xlabel('Elevation angle \theta [deg]');
ylabel('Path loss [dB]');
title('Free-Space Path Loss with Rain Attenuation');
legend('No Attenuation', 'Rain Attenuation: 1%', 'Rain Attenuation: 0.1%', 'Rain Attenuation: 0.01%', 'Rain Attenuation: 0.001%', 'Location', 'best');
grid on;

%% Received Power (Friis Formula, eq 3.8)
P_rxW = EIRPW * G_rx * (lambda / (4*pi) )^2 * (1./d).^2;
P_rx_dBm = 10*log10( P_rxW / (1e-3) ) - Attenuation_list_L_r;
P_rx_dBW = P_rx_dBm - 30;

figure();
plot(theta*180/pi, P_rx_dBm');
xlabel('Elevation angle \theta [deg]');
ylabel('Received power [dBm]');
title('Received Power via Friis Formula');
legend('No Attenuation', 'Rain Attenuation: 1%', 'Rain Attenuation: 0.1%', 'Rain Attenuation: 0.01%', 'Rain Attenuation: 0.001%', 'Location', 'best');
grid on;

%% Signal-to-Noise Ratio (eq 3.91)
SNR_db = P_rx_dBW - F_db - 10*log10(kb*T_receiver*BW);

figure();
plot(theta*180/pi, SNR_db');
xlabel('Elevation angle \theta [deg]');
ylabel('SNR [dB]');
title('Downlink SNR vs. Elevation Angle');
legend('No Attenuation', 'Rain Attenuation: 1%', 'Rain Attenuation: 0.1%', 'Rain Attenuation: 0.01%', 'Rain Attenuation: 0.001%', 'Location', 'best');
grid on;

%% Shannon Channel Capacity
SNR_lin = 10.^(SNR_db/10);
C = BW * log2(1+SNR_lin); % [bit/s]

figure();
plot(theta*180/pi, C');
xlabel('Elevation angle \theta [deg]');
ylabel('Capacity [bit/s]');
title('Shannon Channel Capacity vs. Elevation Angle');
legend('No Attenuation', 'Rain Attenuation: 1%', 'Rain Attenuation: 0.1%', 'Rain Attenuation: 0.01%', 'Rain Attenuation: 0.001%', 'Location', 'best');
grid on;

alpha = asin( (d/(R_E + h_s)) .* cos(theta) );
delta_alpha = abs(alpha(1:end-1) - alpha(2:end));
delta_tn = T_rev * (delta_alpha/(2*pi));

transmitted_data = sum(((C(:,1:end-1) + C(:,2:end))/2) .* delta_tn, 2);

%% Bonus 1.2: Constant Bitrate Optimization
transmission_duration = zeros(1, 91);
transmission_duration(1, 1) = sum(delta_tn);
for i = 1:90
    transmission_duration(1, i+1) = transmission_duration(i) - 2*delta_tn(i);
end
transmitted_data_test(:, 1:91) = C(:, 1:91) .* transmission_duration;
[~, AngleMaxDataTransmitted] = max(transmitted_data_test, [], 2);

figure();
plot(0:90, transmitted_data_test');
xlabel('Minimum cutoff elevation angle \theta [deg]');
ylabel('Transmitted data per pass [bits]');
title('Constant-Bitrate Pass Volume vs. Cutoff Elevation Angle');
legend('No Attenuation', 'Rain Attenuation: 1%', 'Rain Attenuation: 0.1%', 'Rain Attenuation: 0.01%', 'Rain Attenuation: 0.001%', 'Location', 'best');
grid on;

total_data = zeros(5, 1);
for i = 1:length(AngleMaxDataTransmitted)
    total_data(i) = transmitted_data_test(i, AngleMaxDataTransmitted(i));
end

fprintf('Step 1 Summary (No rain, 1%%, 0.1%%, 0.01%%, 0.001%%):\n');
fprintf('  Variable-rate pass volume [Gbit]: %s\n', sprintf('%.4f ', transmitted_data/1e9));
fprintf('  Constant-rate peak volume [Gbit]: %s (optimal cutoff deg: %s)\n', ...
    sprintf('%.4f ', total_data/1e9), sprintf('%d ', AngleMaxDataTransmitted - 1));

%% Export Figures (450 DPI PNG + Vector PDF)
figRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'docs', 'figures');
pngDir = fullfile(figRoot, 'png');
pdfDir = fullfile(figRoot, 'pdf');
if ~exist(pngDir, 'dir'), mkdir(pngDir); end
if ~exist(pdfDir, 'dir'), mkdir(pdfDir); end
figs = findall(0, 'Type', 'figure');
if ~isempty(figs)
    [~, ord] = sort([figs.Number]);
    figs = figs(ord);
    for fi = 1:numel(figs)
        exportgraphics(figs(fi), fullfile(pngDir, sprintf('step1_fig%02d.png', fi)), 'Resolution', 450, 'BackgroundColor', 'white');
        exportgraphics(figs(fi), fullfile(pdfDir, sprintf('step1_fig%02d.pdf', fi)), 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end
