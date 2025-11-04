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

% Initialisation des matrices pour stocker les résultats
Erreur_relative_1 = zeros(1, length(N_values));
Erreur_relative_2 = zeros(1, length(N_values));
Borne_Sup = zeros(1, length(N_values));
Borne_Sup2 = zeros(1, length(N_values));
S_value1=zeros(1, length(N_values));
S_value2=zeros(1, length(N_values));



%% 

for j = 1:length(N_values)
        S_value = S(x_values(1), N_values(j));
        Erreur_relative_1(j) = abs(exp(x_values(1)) - S_value) / exp(x_values(1));
        Borne_Sup(1, j) = exp(-x_values(1))*abs(x_values(1))^(N_values(j) + 1) / factorial(N_values(j) + 1);
        S_value1(j)=S_value;
end

% disp(S_value1)


for j = 1:length(N_values)
        S_value = S(x_values(2), N_values(j));
        Erreur_relative_2(j) = abs(exp(x_values(2)) - S_value) / exp(x_values(2));
        Borne_Sup2(1, j) = abs(x_values(2))^(N_values(j) + 1) / factorial(N_values(j) + 1);  
        S_value2(j)=S_value;
end
% disp(S_value2)

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
S_single = @(x, N) single(sum(single(x).^(0:N) ./ single(factorial(0:N))));

% Calcul des sommes partielles en précision double pour référence
% x_val = -5;
% x_val = 5;
N_val = 0:50;

Sabs_double = zeros(1, length(N_val));
Sabs_double2 = zeros(1, length(N_val));
S_double = zeros(1, length(N_val));
S_double2 = zeros(1, length(N_val));

for j = 1:length(N_val)
    S_double(j) = S_single(x_values(1), N_val(j));
    Sabs_double(j)= abs(S_double(j));
end

for j = 1:length(N_val)
    S_double2(j) = S_single(x_values(2), N_val(j));
    Sabs_double2(j)= abs(S_double2(j));
end


% Calcul des sommes partielles en ordre croissant et décroissant

Sc_val = sort(Sabs_double);
Sd_val = sort(Sabs_double,'descend');

Sc_val2 = sort(Sabs_double2);
Sd_val2 = sort(Sabs_double2,'descend');


% Affichage des graphiques semilogy pour x=-5
figure;
semilogy(N_val, abs((Sc_val - S_value1))./abs(S_value1), 'b', 'LineWidth', 2);
hold on;
semilogy(N_val, abs((Sd_val - S_value1))./abs(S_value1), '--r', 'LineWidth', 2);

xlabel('N');
ylabel('Erreur absolue');
legend('Sc avec x=-5', 'Sd avec x=-5');
title('Erreur absolue de Sc(x;N) et Sd(x;N) par rapport à S(x;N)');



% Affichage des graphiques semilogy pour x=5
figure;
semilogy(N_val, abs((Sc_val2 - S_value2))./abs(S_value2), 'b', 'LineWidth', 2);
hold on;
semilogy(N_val, abs((Sd_val2 - S_value2))./abs(S_value2), '--r', 'LineWidth', 2);

xlabel('N');
ylabel('Erreur absolue');
legend('Sc avec x=5', 'Sd avec x=5');
title('Erreur absolue de Sc(x;N) et Sd(x;N) par rapport à S(x;N)');










