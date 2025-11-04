clear all;
clc
%%  Fonction pour calculer S(x;N) en précision simple
S_single = @(x, N) sum(single(x).^(0:N) ./ single(factorial(0:N)));

% Calcul des sommes partielles en précision double pour référence
x_val = -5;
N_val = 0:50;

S_double = zeros(1, length(N_val));
for j = 1:length(N_val)
    S_double(j) = S_single(x_val, N_val(j));
end

% Calcul des sommes partielles en ordre croissant et décroissant
Sc_val = zeros(length(N_val), length(N_val));
Sd_val = zeros(length(N_val), length(N_val));

for i = 1:length(N_val)
    for j = 1:N_val(i)+1
        Sc_val(i, j) = S_single(x_val, j-1);
        Sd_val(i, j) = S_single(x_val, N_val(i)+1-j);
    end
end

% Affichage des graphiques semilogy
figure;
semilogy(N_val, abs(Sc_val - S_double), '-o', 'LineWidth', 2);
hold on;
semilogy(N_val, abs(Sd_val - S_double), '-^', 'LineWidth', 2);

xlabel('N');
ylabel('Erreur absolue');
legend('Sc(x;N) - S(x;N)', 'Sd(x;N) - S(x;N)');
title('Erreur absolue de Sc(x;N) et Sd(x;N) par rapport à S(x;N)');