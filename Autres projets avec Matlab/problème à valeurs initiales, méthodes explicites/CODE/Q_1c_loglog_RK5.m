%% TP2 MAT6470_2024
clear all; 
clc;
%%

x = [0, -sqrt(3)/4, sqrt(3)/4];
y = [1/2, -1/4, -1/4];

J = 3;
gamma = 4 * pi^2;

% Valeurs de h
alpha = 2.^(1:23);
h_values = 1./alpha(1:23);

% Initialisation des vecteurs pour stocker les résultats
errors_inf_rk5 = zeros(size(h_values));
errors_2_rk5 = zeros(size(h_values));

% Calcul des erreurs pour différentes valeurs de h
for k = 1:length(h_values)
    h = h_values(k);
    timesteps = round(1 / h) + 1;

    % Calcul avec la méthode de RK5
    z_val_num_rk5 = RK5_method(@Res_Syst, [transpose(x); transpose(y)], h, timesteps, J, gamma);

    % Calcul des erreurs d'ordre 2 et infini par rapport à la solution exacte à t=1
    exact_solution = z_val_num_rk5(:, 1);
    
    % Errors for RK5
    errors_2_rk5(k) = sqrt(sum((z_val_num_rk5(:, end) - exact_solution).^2)); % Erreur d'ordre 2
    errors_inf_rk5(k) = norm(z_val_num_rk5(:, end) - exact_solution, inf); % Erreur infinie
end

% Affichage du graphique log-log pour RK5
figure;
loglog(h_values, errors_inf_rk5, 'b-*', 'DisplayName', 'Erreur infinie (RK5)');
hold on;
loglog(h_values, errors_2_rk5, 'r-.', 'DisplayName', 'Erreur d''ordre 2 (RK5)');
loglog(h_values, h_values.^5, 'k--.', 'DisplayName', 'h^5');
xlabel('h');
ylabel('Erreur');
title('Convergence d''ordre 2 et d''ordre infini avec la méthode de RK5');
legend('Location', 'Best');
grid on;
hold off;

% Utilisation de polyfit pour déterminer la pente pour RK5
order_rk5_2 = polyfit(log(h_values(6:13)), log(errors_2_rk5(6:13)), 1);

disp('Ordre de convergence pour l''erreur d''ordre 2 (RK5):');
disp(order_rk5_2(1));

%% Système d'équations
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

% Runge-Kutta 5 (RK5) method
function z_val_rk5 = RK5_method(ode_func, z0, h, timesteps, J, gamma)
    z_val_rk5 = zeros(J + J, timesteps);
    z_val_rk5(:, 1) = z0;

    for i = 2:timesteps
        k1 = h * Res_Syst(z_val_rk5(:, i - 1), J, gamma);
        k2 = h * Res_Syst(z_val_rk5(:, i - 1) + k1 / 4, J, gamma);
        k3 = h * Res_Syst(z_val_rk5(:, i - 1) + 3 * k1 / 32 + 9 * k2 / 32, J, gamma);
        k4 = h * Res_Syst(z_val_rk5(:, i - 1) + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197, J, gamma);
        k5 = h * Res_Syst(z_val_rk5(:, i - 1) + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104, J, gamma);
        k6 = h * Res_Syst(z_val_rk5(:, i - 1) - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40, J, gamma);

        z_val_rk5(:, i) = z_val_rk5(:, i - 1) + 16 * k1 / 135 + 6656 * k3 / 12825 + 28561 * k4 / 56430 - 9 * k5 / 50 + 2 * k6 / 55;
end


end