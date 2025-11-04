%% Convergence avec Heun
clear all;
clc
%% 
x = [0, -sqrt(3)/4, sqrt(3)/4];
y = [1/2, -1/4, -1/4];

J = 3;
gamma = 4 * pi^2;

% Valeurs de h
alpha = 2.^(1:23);
h_values = 1./alpha(1:23);

% Initialisation des vecteurs pour stocker les résultats
errors_inf = zeros(size(h_values));
errors_2 = zeros(size(h_values));

% Calcul des erreurs pour différentes valeurs de h
for k = 1:length(h_values)
    h = h_values(k);
    timesteps = round(1 / h) + 1;

    % Calcul avec la méthode de Heun
    z_val_num = zeros(J + J, timesteps);
    z_val_num(:, 1) = [transpose(x); transpose(y)];

    for i = 2:timesteps
        % Predictor Step (Euler Method)
        pred = z_val_num(:, i - 1) + h * Res_Syst(z_val_num(:, i - 1), J, gamma);

        % Corrector Step
        z_val_num(:, i) = z_val_num(:, i - 1) + 0.5 * h * (Res_Syst(z_val_num(:, i - 1), J, gamma) + Res_Syst(pred, J, gamma));
    end

    % Calcul des erreurs d'ordre 2 et infini par rapport à la solution exacte à t=1
    exact_solution = z_val_num(:, 1);
    errors_2(k) = sqrt(sum((z_val_num(:, end) - exact_solution).^2)); % Erreur d'ordre 2
    errors_inf(k) = norm(z_val_num(:, end) - exact_solution, inf); % Erreur infinie
end

% Affichage du graphique log-log
figure;
loglog(h_values, errors_inf, 'b-*', 'DisplayName', 'Erreur infinie');
hold on;
loglog(h_values, errors_2, 'r-.', 'DisplayName', 'Erreur d''ordre 2');
loglog(h_values, h_values.*h_values, 'k--.', 'DisplayName', 'h^2');
xlabel('h');
ylabel('Erreur');
title('Ereeur avec la méthode de Heun');
legend('Location', 'Best');
grid on;
hold off;

% Utilisation de polyfit pour déterminer la pente
order_RK2 = polyfit(log(h_values(18:20)), log(errors_2(18:20)), 1);

disp('Ordre de convergence pour l''erreur d''ordre 2:');
disp(order_RK2(1));


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