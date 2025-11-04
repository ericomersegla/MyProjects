clear all;
clc
%% 
% Fonction pour calculer S(x;N)
S = @(x, N) sum(x.^(0:N) ./ factorial(0:N));
%% 




% Valeurs de x
x_values = [-5, 5];

% Valeurs de N
N_values = 0:50; 


    for j = 1:length(N_values)
        S_value_1 = S(-5, N_values(j));
        disp(S_value_1)
    end

disp('kjkjjjhjjjjj')

    for j = 1:length(N_values)
        S_value_2 = S(5, N_values(j));
        disp(S_value_2)
    end


% Initialisation des matrices pour stocker les résultats
Erreur_relative = zeros(length(x_values), length(N_values));
Borne_Sup = zeros(length(x_values), length(N_values));
%% 
% Calcul et stockage des erreurs relatives et de la borne superieure
for i = 1:length(x_values)
    for j = 1:length(N_values)
        S_value = S(x_values(i), N_values(j));
        Erreur_relative(i, j) = abs(exp(x_values(i)) - S_value) / exp(x_values(i));
        Borne_Sup(i, j) = abs(x_values(i))^(N_values(j) + 1) / factorial(N_values(j) + 1);
    end
end
%% 
% Tracé du graphique semilogy
figure;
semilogy(N_values, Erreur_relative(1, :), 'b', 'LineWidth', 2); 
hold on;
semilogy(N_values, Erreur_relative(2, :), 'r', 'LineWidth', 2);
semilogy(N_values, Borne_Sup(1, :), '--b', 'LineWidth', 2);
semilogy(N_values, Borne_Sup(2, :), '--r', 'LineWidth', 2);
xlabel('N');
ylabel('Erreur relative');
legend('x = -5', 'x = 5', 'Borne Superieure théorique pour x = -5', 'Borne Superieure théorique pour x = 5');
title('Erreur relative de S(x;N) par rapport à exp(x)');
