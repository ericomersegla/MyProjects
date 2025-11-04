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
errors_inf_rk4 = zeros(size(h_values));
errors_2_rk4 = zeros(size(h_values));

% Calcul des erreurs pour différentes valeurs de h
for k = 1:length(h_values)
    h = h_values(k);
    timesteps = round(1 / h) + 1;

    % Runge-Kutta 4 (RK4) method
    z_val_num_rk4 = zeros(J + J, timesteps);
    z_val_num_rk4(:, 1) = [transpose(x); transpose(y)];

    for i = 2:timesteps
        k1 = h * Res_Syst(z_val_num_rk4(:, i - 1), J, gamma);
        k2 = h * Res_Syst(z_val_num_rk4(:, i - 1) + 0.5 * k1, J, gamma);
        k3 = h * Res_Syst(z_val_num_rk4(:, i - 1) + 0.5 * k2, J, gamma);
        k4 = h * Res_Syst(z_val_num_rk4(:, i - 1) + k3, J, gamma);

        z_val_num_rk4(:, i) = z_val_num_rk4(:, i - 1) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    end

    % Calcul des erreurs d'ordre 2 et infini par rapport à la solution exacte à t=1
    exact_solution = z_val_num_rk4(:, 1);
    
    % Errors for RK4
    errors_2_rk4(k) = sqrt(sum((z_val_num_rk4(:, end) - exact_solution).^2)); % Erreur d'ordre 2
    errors_inf_rk4(k) = norm(z_val_num_rk4(:, end) - exact_solution, inf); % Erreur infinie


end

% Affichage du graphique log-log pour RK4
figure;
loglog(h_values, errors_inf_rk4, 'r', 'DisplayName', 'Erreur infinie (RK4)');
hold on;
loglog(h_values, errors_2_rk4, 'b', 'DisplayName', 'Erreur d''ordre 2 (RK4)');
loglog(h_values, h_values.^4, 'k--', 'DisplayName', 'h^4');
xlabel('h');
ylabel('Erreur');
title('Convergence d''ordre 2 et d''ordre infini avec la méthode de RK4');
legend('Location', 'Best');
grid on;
hold off;



% Utilisation de polyfit pour déterminer la pente
order_RK4 = polyfit(log(h_values(11:14)), log(errors_2_rk4(11:14)), 1);

disp('Ordre de convergence pour l''erreur d''ordre 2 de RK4:');
disp(order_RK4(1));


% Function for Runge-Kutta 4 method
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