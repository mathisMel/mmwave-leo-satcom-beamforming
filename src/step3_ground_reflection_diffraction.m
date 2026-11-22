%% 26 GHz mmWave LEO SatCom — Step 3: Ground Reflection & Building Diffraction (NLOS/LOS)
% Course: ELEC-H415 (Wireless Communications / Propagation Modeling)
%
% Models an urban street-canyon receiver flanked by 10 m tall dielectric walls
% (x = +/-10 m, epsilon_r = 4) over theta = 0..180 deg:
%   - theta in [0, 19) and (162, 180] deg: single-edge Fresnel diffraction (NLOS)
%   - theta in [19, 45] and [135, 161] deg: specular Fresnel reflection off building walls
%   - theta in (45, 135) deg: unobstructed LOS reception
% Evaluates 1-antenna baseline plus 4, 8, and 16-element isotropic & patch ULAs.
%
% Runtime note: Default grid step pas = 0.2 deg runs in ~7-20 s. Setting
% pas = 0.02 deg matches the high-resolution report sweep but is slow due to
% full 9001x9001x3 steering arrays.

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end

%% Parameters
fc = 26e9;                 % Carrier frequency [Hz]
c = 3e8;                   % Speed of light [m/s]
lambda = c/fc;             % Wavelength [m]
beta = 2*pi/ lambda;       % Wavenumber [1/m]
omega = 2*pi*fc;           % Angular frequency [rad/s]
epsilon_r = 4;             % Relative permittivity of building facade
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

pas = 0.2;                 % Angular grid resolution [deg] (use 0.02 for report high-res grid)
theta = (0:pas:180)*pi/180;% Elevation angle grid [rad]
d = sqrt( ((R_E+h_s)^2) - ( (R_E^2) * (cos(theta).^2) ) ) - (R_E.*sin(theta)); % Slant range [m]
idx45 = round(45/pas) + 1; % Index corresponding to 45 deg wall-top ray

% ULA parameters
antenna_space = lambda/2;
U_0 = 1;
size_array = [4; 8; 16];                           % Number of antennas in the array
size_array = permute(size_array, [3 2 1]);
theta_m = [90, 45, 30, 0] * pi/180;
delta = (beta*antenna_space.*cos(theta))';
psi = beta*antenna_space.*cos(theta) - delta;

% Knife-edge diffraction parameters
h = 10;                                            % Wall height [m]
d2 = zeros(1, 92/pas);
nu = zeros(1, 92/pas);
d1 = 10;                                           % Distance from wall to RX [m]
d2(1:46/pas) = d(1:46/pas) .* cos(theta(1:46/pas)) - 10;
d2((46+pas)/pas:end) = flip(d2(1:46/pas));
nu(1:46/pas) = (1-tan(theta(1:46/pas)))*h .* sqrt( (2/lambda) * (1./d2(1:46/pas) + 1/d1) );
nu((46+pas)/pas:end) = flip(nu(1:46/pas));
norm__squared_F_nu = -6.9 - 20*log10( sqrt( (nu-0.1).^2 + 1 ) + nu - 0.1); % Eq 3.70

% Total transmitted data parameters
alpha = asin( (d/(R_E + h_s)) .* cos(theta) );
delta_alpha = alpha(1:end-1) - alpha(2:end);
delta_tn = T_rev * (delta_alpha / (2*pi));

% Street canyon geometry: walls at x = +/-10 m, height 10 m
coef_reflex = ( cos(theta(20/pas:46/pas)) - ( sqrt(epsilon_r) * sqrt( 1 - ( (1/epsilon_r) * sin(theta(20/pas:46/pas)).^2) ) ) ) ./ ...
              ( cos(theta(20/pas:46/pas)) + ( sqrt(epsilon_r) * sqrt( 1 - ( (1/epsilon_r) * sin(theta(20/pas:46/pas)).^2) ) ) );
theta_deg = theta * 180 / pi;

%% Step 1 Baseline: Single Isotropic Antenna (NLOS + LOS)
G_rx = 1;
A_e_rx = ( lambda^2 / (4*pi) ) * G_rx;
E = ( sqrt(60*EIRPW) * exp(-1j*beta*d) ) ./ d;
P_rxW_WithoutReflexion = (A_e_rx * arrayfun(@norm, E).^2 ) / (2* Z_0);

E(20/pas:46/pas) = E(20/pas:46/pas) .* coef_reflex;
E(136/pas: end-((18+pas)/pas)) = E(136/pas: end-((18+pas)/pas)) .* flip(coef_reflex);
P_rxW_With_Reflexion = (A_e_rx * arrayfun(@norm, E).^2 ) / (2* Z_0);

P_rx_dBW_WithoutReflexion = 10*log10( P_rxW_WithoutReflexion );
P_rx_dBW_WithReflexion = 10*log10( P_rxW_With_Reflexion );

Reflex_dbWStart = abs(P_rx_dBW_WithReflexion(1:46/pas) - P_rx_dBW_WithoutReflexion(1:46/pas));
Reflex_dbWEnd = abs(P_rx_dBW_WithReflexion(136/pas:end) - P_rx_dBW_WithoutReflexion(136/pas:end));
indices_No_ReflexStart = find( Reflex_dbWStart == 0);
indices_No_ReflexEnd = find(Reflex_dbWEnd == 0);

P_rx_dBW_WithReflexion(indices_No_ReflexStart) = P_rx_dBW_WithoutReflexion(indices_No_ReflexStart) + norm__squared_F_nu(indices_No_ReflexStart);
P_rx_dBW_WithReflexion(136/pas+indices_No_ReflexEnd-1) = P_rx_dBW_WithoutReflexion(136/pas+indices_No_ReflexEnd-1) + norm__squared_F_nu(46/pas+indices_No_ReflexEnd);
P_rx_dBm_With_diffraction = P_rx_dBW_WithReflexion + 30;

SNR_db = P_rx_dBW_WithReflexion - F_db - 10*log10(k*T_receiver*BW);
SNR_lin = 10.^(SNR_db/10);
C = BW * log2(1+SNR_lin);
transmitted_data_1Antenna_iso = sum(((C(1:end-1) + C(2:end))/2) .* delta_tn, 2);

figure;
subplot(2,1,1);
plot(theta_deg, P_rx_dBm_With_diffraction);
title('Received Power (NLOS + LOS) — Single Isotropic Antenna');
xlabel('Elevation angle \theta [deg]'); ylabel('P_{rx} [dBm]'); grid on;

subplot(2,1,2);
plot(theta_deg, C);
xlabel('Elevation angle \theta [deg]'); ylabel('Capacity [bit/s]');
title('Channel Capacity (NLOS + LOS) — Single Isotropic Antenna'); grid on;

%% Isotropic ULA Step 2 (4, 8, 16 Antennas)
U_theta = U_0;
psi(psi == 0) = 1e-15;
psi(4, length(theta)) = 1e-15;

Array_Factor = ( exp(1j*(size_array-1).*psi/2) .* sin(size_array.*psi/2) ) ./ (sin(psi/2));
U_array = arrayfun(@abs, Array_Factor).^2 .* U_theta;
sin_theta = sin(theta);
P_rx = trapz(theta, U_array .* sin_theta, 2);

[row, ~] = size(U_array);
for i = 1:row
    if i <= 19/pas || i >= 163/pas
        U_array(i, 1, :) = U_array(i, idx45, :);
    elseif (i > 19/pas && i < 47/pas) || (i > 135/pas && i < 163/pas)
        U_array(i, 1, :) = U_array(i, end-i+1, :);
    else
        U_array(i, 1, :) = U_array(i, i, :);
    end
end
U_array = U_array(:, 1, :);
Gain_rx = U_array ./ ( (1/2) * P_rx );
Gain_rx = permute(Gain_rx, [3 1 2]);

A_e_rx = ( lambda^2 / (4*pi) ) * Gain_rx;
E = ( sqrt(60*EIRPW) * exp(-1j*beta*d) ) ./ d;
P_rxW_WithoutReflexion = (A_e_rx .* arrayfun(@norm, E).^2 ) / (2* Z_0);
E(20/pas:46/pas) = E(20/pas:46/pas) .* coef_reflex;
E(136/pas: end-((18+pas)/pas)) = E(136/pas: end-((18+pas)/pas)) .* flip(coef_reflex);
P_rxW_With_Reflexion = (A_e_rx .* arrayfun(@norm, E).^2 ) / (2* Z_0);

P_rx_dBW_WithoutReflexion = 10*log10( P_rxW_WithoutReflexion );
P_rx_dBW_WithReflexion = 10*log10( P_rxW_With_Reflexion );

Reflex_dbWStart = abs(P_rx_dBW_WithReflexion(1,1:46/pas) - P_rx_dBW_WithoutReflexion(1,1:46/pas));
Reflex_dbWEnd = abs(P_rx_dBW_WithReflexion(1,136/pas:end) - P_rx_dBW_WithoutReflexion(1,136/pas:end));
indices_No_ReflexStart = find( Reflex_dbWStart == 0);
indices_No_ReflexEnd = find(Reflex_dbWEnd == 0);

P_rx_dBW_WithReflexion(:,indices_No_ReflexStart) = P_rx_dBW_WithoutReflexion(:,indices_No_ReflexStart) + norm__squared_F_nu(indices_No_ReflexStart);
P_rx_dBW_WithReflexion(:,136/pas+indices_No_ReflexEnd-1) = P_rx_dBW_WithoutReflexion(:,136/pas+indices_No_ReflexEnd-1) + norm__squared_F_nu(47/pas+indices_No_ReflexEnd-1);
P_rx_dBm_With_diffraction = P_rx_dBW_WithReflexion + 30;

SNR_db = P_rx_dBW_WithReflexion - F_db - 10*log10(k*T_receiver*BW);
SNR_lin = 10.^(SNR_db/10);
C = BW * log2(1+SNR_lin);
transmitted_data_iso = sum(((C(:,1:end-1,:) + C(:,2:end,:))/2) .* delta_tn, 2);
transmitted_data_iso = squeeze(transmitted_data_iso);
max1 = max(C, [], 2);

figure;
subplot(2,1,1);
plot(theta_deg, P_rx_dBm_With_diffraction');
title('Received Power (NLOS + LOS) — 4, 8, 16 Isotropic Antennas');
xlabel('Elevation angle \theta [deg]'); ylabel('P_{rx} [dBm]');
legend('4 antennas', '8 antennas', '16 antennas'); grid on;

subplot(2,1,2);
plot(theta_deg, C');
title('Channel Capacity (NLOS + LOS) — 4, 8, 16 Isotropic Antennas');
xlabel('Elevation angle \theta [deg]'); ylabel('Capacity [bit/s]');
legend('4 antennas', '8 antennas', '16 antennas'); grid on;

%% Patch ULA Step 2 (4, 8, 16 Antennas)
U_theta = U_0 * sin(theta).^2;
psi(psi == 0) = 1e-15;
psi(4, length(theta)) = 1e-15;

Array_Factor = ( exp(1j*(size_array-1).*psi/2) .* sin(size_array.*psi/2) ) ./ (sin(psi/2));
U_array = arrayfun(@abs, Array_Factor).^2 .* U_theta;
P_rx = trapz(theta, U_array .* sin(theta), 2);

[row, ~] = size(U_array);
for i = 1:row
    if i <= 19/pas || i >= 163/pas
        U_array(i, 1, :) = U_array(i, idx45, :);
    elseif (i > 19/pas && i < 47/pas) || (i > 135/pas && i < 163/pas)
        U_array(i, 1, :) = U_array(i, end-i+1, :);
    else
        U_array(i, 1, :) = U_array(i, i, :);
    end
end
U_array = U_array(:, 1, :);
Gain_rx = U_array ./ ( (1/2) * P_rx );
Gain_rx = permute(Gain_rx, [3 1 2]);

A_e_rx = ( lambda^2 / (4*pi) ) * Gain_rx;
E = ( sqrt(60*EIRPW) * exp(-1j*beta*d) ) ./ d;
P_rxW_WithoutReflexion = (A_e_rx .* arrayfun(@norm, E).^2 ) / (2* Z_0);
E(20/pas:46/pas) = E(20/pas:46/pas) .* coef_reflex;
E(136/pas: end-((18+pas)/pas)) = E(136/pas: end-((18+pas)/pas)) .* flip(coef_reflex);
P_rxW_With_Reflexion = (A_e_rx .* arrayfun(@norm, E).^2 ) / (2* Z_0);
P_rxW_With_Reflexion(:, length(theta), :) = 0;
P_rxW_WithoutReflexion(:, length(theta), :) = 0;

P_rx_dBW_WithoutReflexion = 10*log10( P_rxW_WithoutReflexion );
P_rx_dBW_WithReflexion = 10*log10( P_rxW_With_Reflexion );

P_rx_dBW_WithReflexion(:,indices_No_ReflexStart) = P_rx_dBW_WithoutReflexion(:,indices_No_ReflexStart) + norm__squared_F_nu(indices_No_ReflexStart);
P_rx_dBW_WithReflexion(:,136/pas+indices_No_ReflexEnd-1) = P_rx_dBW_WithoutReflexion(:,136/pas+indices_No_ReflexEnd-1) + norm__squared_F_nu(47/pas+indices_No_ReflexEnd-1);
P_rx_dBm_With_diffraction = P_rx_dBW_WithReflexion + 30;

SNR_db = P_rx_dBW_WithReflexion - F_db - 10*log10(k*T_receiver*BW);
SNR_lin = 10.^(SNR_db/10);
C = BW * log2(1+SNR_lin);
transmitted_data_patch = sum(((C(:,1:end-1,:) + C(:,2:end,:))/2) .* delta_tn, 2);
transmitted_data_patch = squeeze(transmitted_data_patch);
max2 = max(C, [], 2);

figure;
subplot(2,1,1);
plot(theta_deg, P_rx_dBm_With_diffraction');
title('Received Power (NLOS + LOS) — 4, 8, 16 Patch Antennas');
xlabel('Elevation angle \theta [deg]'); ylabel('P_{rx} [dBm]');
legend('4 antennas', '8 antennas', '16 antennas'); grid on;

subplot(2,1,2);
plot(theta_deg, C');
title('Channel Capacity (NLOS + LOS) — 4, 8, 16 Patch Antennas');
xlabel('Elevation angle \theta [deg]'); ylabel('Capacity [bit/s]');
legend('4 antennas', '8 antennas', '16 antennas'); grid on;

fprintf('Step 3 Summary (pas = %.2f deg):\n', pas);
fprintf('  1-antenna Iso pass data   : %.4f Gbit\n', transmitted_data_1Antenna_iso/1e9);
fprintf('  [4, 8, 16] Iso pass data  : %s Gbit | Peak C: %s Mbit/s\n', ...
    sprintf('%.3f ', transmitted_data_iso/1e9), sprintf('%.2f ', max1/1e6));
fprintf('  [4, 8, 16] Patch pass data: %s Gbit | Peak C: %s Mbit/s\n', ...
    sprintf('%.3f ', transmitted_data_patch/1e9), sprintf('%.2f ', max2/1e6));

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
        exportgraphics(figs(fi), fullfile(pngDir, sprintf('step3_fig%02d.png', fi)), 'Resolution', 450, 'BackgroundColor', 'white');
        exportgraphics(figs(fi), fullfile(pdfDir, sprintf('step3_fig%02d.pdf', fi)), 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end
