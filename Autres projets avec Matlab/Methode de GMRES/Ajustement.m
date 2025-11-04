
n_values = 2:200; % Valeurs de n à tester
cpu_time = zeros(length(n_values),1); % Pour stocker les temps CPU moyens

for j = 1:length(n_values)
    n = n_values(j);
    [~, ~, time] = GMRES_Method(A, b, tol,n);
    cpu_time (j)=sum(time)/199;
end


plot(n_values,cpu_time,'r-*');
xlabel('n');
ylabel('Temps CPU moyen');
title('Temps CPU moyen pour la factorisation QR en fonction de n');
grid on
% 

% Ajustement linéaire
coefficients = polyfit(n_values, cpu_time, 1);
slope = coefficients(1); % Pente de la droite ajustée
intercept = coefficients(2); % Ordonnée à l'origine de la droite ajustée

% Création de la droite ajustée
fit_line = slope * n_values + intercept;

% Ajout de la droite ajustée au graphique
hold on;
plot(n_values, fit_line, 'b--', 'LineWidth', 1);
legend('Données', 'Ajustement linéaire', 'Location', 'northwest');
hold off;
