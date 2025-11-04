clear all;
clc
%% 
%2- Fonction pour calculer S(x;N)
S = @(x, N) sum(x.^(0:N) ./ factorial(0:N));
%%
% Valeurs de x
x_values = [-5, 5];

% Valeurs de N
N_values = 0:50;

% Calcul et stockage des erreurs relatives et de la borne superieure
% for i = 1:length(x_values)
%     for j = 1:length(N_values)
%         S_value = S(x_values(i), N_values(j));
%         Erreur_relative(i, j) = abs(exp(x_values(i)) - S_value) / exp(x_values(i));
%         Borne_Sup(i, j) = abs(x_values(i))^(N_values(j) + 1) / factorial(N_values(j) + 1);
%         disp(S_value)
%     end
% end

% Initialisation des matrices pour stocker les résultats
Erreur_relative_1 = zeros(1, length(N_values));
Erreur_relative_2 = zeros(1, length(N_values));
Borne_Sup = zeros(1, length(N_values));
Borne_Sup2 = zeros(1, length(N_values));
%% 
% 
% for i = 1:length(x_values)
%     for j = 1:length(N_values)
%         S_value = S(x_values(i), N_values(j));
%         Erreur_relative(i, j) = abs(exp(x_values(i)) - S_value) / exp(x_values(i));
%         Borne_Sup(1, j) = abs(x_values(1))^(N_values(j) + 1) / factorial(N_values(j) + 1);
%         Borne_Sup2(1, j) = abs(x_values(2))^(N_values(j) + 1) / factorial(N_values(j) + 1);
%         disp(S_value)
%     end
% end

for j = 1:length(N_values)
        S_value = S(x_values(1), N_values(j));
        Erreur_relative_1(j) = abs(exp(x_values(1)) - S_value) / exp(x_values(1));
        Borne_Sup(1, j) = abs(x_values(1))^(N_values(j) + 1) / factorial(N_values(j) + 1);
        disp(S_value)
    end


for j = 1:length(N_values)
        S_value = S(x_values(2), N_values(j));
        Erreur_relative_2(j) = abs(exp(x_values(2)) - S_value) / exp(x_values(2));
        Borne_Sup2(1, j) = exp(-x_values(2)).abs(x_values(1))^(N_values(j) + 1) / factorial(N_values(j) + 1);
        disp(S_value)
    end




%% 
% Tracé du graphique semilogy
figure;
semilogy(N_values, Erreur_relative_1, 'b', 'LineWidth', 2); 
hold on;
semilogy(N_values, Erreur_relative_2, 'r', 'LineWidth', 2);
semilogy(N_values, Borne_Sup, '--b', 'LineWidth', 2);
semilogy(N_values, Borne_Sup2, '--r', 'LineWidth', 2);
xlabel('N');
ylabel('Erreur relative');
legend('x = -5', 'x = 5', 'Borne Superieure théorique pour x = -5', 'Borne Superieure théorique pour x = 5');
title('Erreur relative de S(x;N) par rapport à exp(x)');
%%  3-Fonction pour calculer S(x;N) en précision simple
S_single = @(x, N) sum(single(x).^(0:N) ./ single(factorial(0:N)));

% Calcul des sommes partielles en précision double pour référence
% x_val = -5;
x_val = 5;
N_val = 0:50;

Sabs_double = zeros(1, length(N_val));
S_double = zeros(1, length(N_val));
for j = 1:length(N_val)
    S_double(j) = S_single(x_val, N_val(j));
    Sabs_double(j)= abs(S_double(j));
end

% Calcul des sommes partielles en ordre croissant et décroissant
% Sc_val = zeros(length(N_val), length(N_val));
% Sd_val = zeros(length(N_val), length(N_val));

Sc_val = sort(Sabs_double);
Sd_val = sort(Sabs_double,'descend');


for i = 1:length(N_val)
    for j = 1:N_val(i)+1
        Sc_val(i, j) = S_single(x_val, j-1);
        Sd_val(i, j) = S_single(x_val, N_val(i)+1-j);
    end
end
A=abs(Sc_val - S_double)/S_double;
B=abs(Sd_val - S_double)/S_double;

% Affichage des graphiques semilogy
figure;
semilogy(N_val, A, '--o', 'LineWidth', 2);
hold on;
semilogy(N_val, B, '--^', 'LineWidth', 2);

xlabel('N');
ylabel('Erreur absolue');
legend('Sc(x;N) - S(x;N)', 'Sd(x;N) - S(x;N)');
title('Erreur absolue de Sc(x;N) et Sd(x;N) par rapport à S(x;N)');
