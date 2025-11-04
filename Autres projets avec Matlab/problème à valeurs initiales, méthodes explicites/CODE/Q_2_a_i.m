%% TP2 MAT6470_2024
clear all; 
clc;
%% 
theta = (4 * pi) / 9;
t_0 = 0;
t_f = 9.06;
x = [(3 + sqrt(3) * cos(theta)) / 6, (-3 + sqrt(3) * cos(theta)) / 6, (2 * sqrt(3) / 3) * cos(theta)];
y = [(sqrt(3) / 6) * sin(theta), (sqrt(3) / 6) * sin(theta), (2 * sqrt(3) / 3) * sin(theta)];

J = 3;

% Choisir trois valeurs différentes pour gamma
gamma_values = [2, 1, 3];

figure;
hold on

for k = 1:length(gamma_values)
    gamma = gamma_values(k);

    h_max = 0.2;
    h = h_max;
    timesteps = ceil((t_f - t_0) / h) + 1; % Correction pour s'assurer que le nombre de pas de temps est un entier
    z_val_rkf = zeros(J + J, timesteps);
    z_val_rkf(:, 1) = [transpose(x); transpose(y)];  % Les premières valeurs sont les valeurs initiales

    for i = 2:timesteps
        % Runge-Kutta-Fehlberg (RKF)
        k1 = h * Res_Syst(z_val_rkf(:, i - 1), J, gamma);
        k2 = h * Res_Syst(z_val_rkf(:, i - 1) + k1 / 4, J, gamma);
        k3 = h * Res_Syst(z_val_rkf(:, i - 1) + 3 * k1 / 32 + 9 * k2 / 32, J, gamma);
        k4 = h * Res_Syst(z_val_rkf(:, i - 1) + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197, J, gamma);
        k5 = h * Res_Syst(z_val_rkf(:, i - 1) + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104, J, gamma);
        k6 = h * Res_Syst(z_val_rkf(:, i - 1) - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40, J, gamma);

        % solutions
        Err_rkf = k1 / 360 - 128 * k3 / 4275 - 2197 * k4 / 75240 + k5 / 50 + 2 * k6 / 55;
        z_rkf = z_val_rkf(:, i - 1) + Err_rkf;

        % Calcul de l'erreur
        error = norm(Err_rkf, 'inf');

        % Adaptation du pas de temps
        tolerance = 1e-13 + norm(z_val_rkf(:, i - 1), 'inf') * 5e-10;
        h_new = h * min(max(0.84 * (tolerance / error)^0.25, 0.1), 4);

        % Mise à jour du pas de temps
        h = min(h_new, h_max);

        % Mise à jour des valeurs avec le pas de temps adaptatif
        z_val_rkf(:, i) = z_rkf;
    end

    % Affichage des trajectoires avec RKF et pas de temps adaptatif
    plot3(z_val_rkf(1, :), z_val_rkf(4, :), linspace(t_0, t_f, timesteps), "r-*");
    hold on
    plot3(z_val_rkf(2, :), z_val_rkf(5, :), linspace(t_0, t_f, timesteps), "b-*");
    hold on
    plot3(z_val_rkf(3, :), z_val_rkf(6, :), linspace(t_0, t_f, timesteps), "k-*");
end

xlabel("xi")
ylabel("Yi")
zlabel("Time")
legend('Tourbillon 1 (\gamma = 2)', 'Tourbillon 2 (\gamma = 1)', 'Tourbillon 3 (\gamma = 3)')
title('Trajectoires avec RKF et pas de temps adaptatif')

function dzdt = Res_Syst(z, J, G)
    dzdt = zeros(J + J, 1);

    for i = 1:J
        dxdt = 0;
        dydt = 0;

        for j = 1:J
            if i ~= j
                dxdt = dxdt + G / (2 * pi) * (z(J+i) - z(J+j)) / ((z(i) - z(j))^2 + (z(J+i) - z(J+j))^2);
                dydt = dydt + G / (2 * pi) * (z(i) - z(j)) / ((z(i) - z(j))^2 + (z(J+i) - z(J+j))^2);
            end
        end

        dzdt(i) = -dxdt;
        dzdt(J+i) = dydt;
    end
end