%% Refaire les mêmes calculs en utilisant la méthode RK4.
%% TP2 MAT6470_2024
clear all; 
clc;
%%

x = [0, -sqrt(3)/4, sqrt(3)/4];
y = [1/2, -1/4, -1/4];

J = 3;
gamma = 4 * pi^2;
h = 0.01;
timesteps = 1 / h + 1;
z_val_rk4 = zeros(J + J, timesteps);
z_val_rk4(:, 1) = [transpose(x); transpose(y)];  % Les premières valeurs sont les valeurs initiales

% Méthode de Runge-Kutta d'ordre 4
for i = 2:timesteps
    k1 = h * Res_Syst(z_val_rk4(:, i - 1), J, gamma);
    k2 = h * Res_Syst(z_val_rk4(:, i - 1) + 0.5 * k1, J, gamma);
    k3 = h * Res_Syst(z_val_rk4(:, i - 1) + 0.5 * k2, J, gamma);
    k4 = h * Res_Syst(z_val_rk4(:, i - 1) + k3, J, gamma);

    z_val_rk4(:, i) = z_val_rk4(:, i - 1) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end

% Affichage des trajectoires avec RK4
figure;
plot(z_val_rk4(1,:), z_val_rk4(4,:), "r--");
hold on
plot(z_val_rk4(2,:), z_val_rk4(5,:), "b--");
hold on
plot(z_val_rk4(3,:), z_val_rk4(6,:), "k--");
xlabel("xi")
ylabel("Yi")
legend('Tourbillon 1', 'Tourbillon 2', 'Tourbillon 3')
title('Trajectoires avec la méthode de RK4')
hold off

% Système d'équation pour RK4
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