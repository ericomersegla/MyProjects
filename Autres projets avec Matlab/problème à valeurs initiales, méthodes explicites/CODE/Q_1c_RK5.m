%% Refaire les mêmes calculs en utilisant la méthode RK5.
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
z_val_rk5 = zeros(J + J, timesteps);
z_val_rk5(:, 1) = [transpose(x); transpose(y)];  % Les premières valeurs sont les valeurs initiales

% Méthode de Runge-Kutta d'ordre 5
for i = 2:timesteps
    k1 = h * Res_Syst(z_val_rk5(:, i - 1), J, gamma);
    k2 = h * Res_Syst(z_val_rk5(:, i - 1) + k1 / 4, J, gamma);
    k3 = h * Res_Syst(z_val_rk5(:, i - 1) + 3 * k1 / 32 + 9 * k2 / 32, J, gamma);
    k4 = h * Res_Syst(z_val_rk5(:, i - 1) + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197, J, gamma);
    k5 = h * Res_Syst(z_val_rk5(:, i - 1) + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104, J, gamma);
    k6 = h * Res_Syst(z_val_rk5(:, i - 1) - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40, J, gamma);

    z_val_rk5(:, i) = z_val_rk5(:, i - 1) + 16 * k1 / 135 + 6656 * k3 / 12825 + 28561 * k4 / 56430 - 9 * k5 / 50 + 2 * k6 / 55;
end

% Affichage des trajectoires avec RK5
figure;
plot(z_val_rk5(1,:), z_val_rk5(4,:), "r-*");
hold on
plot(z_val_rk5(2,:), z_val_rk5(5,:), "b-*");
hold on
plot(z_val_rk5(3,:), z_val_rk5(6,:), "k-*");
xlabel("xi")
ylabel("Yi")
legend('Tourbillon 1', 'Tourbillon 2', 'Tourbillon 3')
title('Trajectoires avec RK5')
hold off

% Système d'équation pour RK5
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