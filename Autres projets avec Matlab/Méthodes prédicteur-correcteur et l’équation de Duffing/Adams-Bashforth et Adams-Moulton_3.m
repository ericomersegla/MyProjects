%% Définition des paramètres
clear all
clc

%% Définition des paramètres
g_m = -1;
l = 1;
gamma = 0.5;
A = 0.5;
omega = 1;
% Conditions initiales
alpha_0 = pi/2;
alpha_1 = 0;
n = 2^10;
h = 100./n;
t = 0:h:100;

y0 = [alpha_0; alpha_1]; % Correction des conditions initiales

% Initialisation de y pour stocker les résultats
y = zeros(2,n+1); % Correction de la taille de y

for i = 1:n+1 % Correction de la boucle
    % Appel de la fonction Mafun avec les paramètres corrects
    [t_temp, y_temp] = Cauchy_AB3_AM4(@(t, y) Mafun(t,gamma, y, g_m, A), [0 100],y0, n);  
    
end

figure;
plot(t, y_temp(1,:), 'b-', 'DisplayName', ' theta 1');
% plot(h, Error, 'b-*', 'DisplayName', 'Erreur du predicteur AB3-AM4');
xlabel('temps');
ylabel('theta 1');
title('Construction de theta 1 en fonction du temps');
grid on;
hold off;

figure;
plot(y_temp(1,:), y_temp(2,:), 'r-', 'DisplayName', ' theta 2');
% plot(h, Error, 'b-*', 'DisplayName', 'Erreur du predicteur AB3-AM4');
xlabel('theta 1');
ylabel('theta 2');
title('Construction de theta 2 en fonction du theta 1');
grid on;
hold off;


function derive = Mafun(t,gamma, y, g, A) % Correction des paramètres
    l=1;
    omega = 1;
    derive = [y(2); -gamma * y(2) - (g /l)* (y(1)-((y(1))^3/6))+ A*cos(omega*t) ]; % Ajout du terme A*sin(omega*t) 
end


function [t, y_rk]=rk4(f,tspan, y0, n)
% Paramètres de la méthode
    h=(tspan(2)-tspan(1))/n;
    t = tspan(1):h:tspan(2);
    y_rk= zeros(length(y0),3);
    y_rk( :,1) = y0;
% Utiliser RK4 pour obtenir les valeurs initiales y1 et y2
    for i = 1:2
        k1 = h * feval(f, t(i), y_rk(:,i));
        k2 = h * feval(f, t(i) + 0.5* h , y_rk(:,i) + 0.5*k1);
        k3 = h * feval(f, t(i) + 0.5*h , y_rk(:,i) + 0.5*k2);
        k4 = h * feval(f, t(i) + h, y_rk(:,i) + k3);

        y_rk(:,i+1) = y_rk(:,i) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    end

end



function [t, y] = Cauchy_AB3_AM4(f,tspan, y0, n)
    % Résout numériquement l'équation différentielle avec la méthode AB3-AM4
    
    [t, y_rk]=rk4(f,tspan, y0, n);
    % Paramètres de la méthode
    h=(tspan(2)-tspan(1))/n;
    t = tspan(1):h:tspan(2);
   
    y = zeros(length(y0),n+1);
    % y_pred = zeros(length(y0),n+1);
    y(:,1) = y0;
    y(:,2) = y_rk(:,2);
    y(:,3)= y_rk(:,3);
    
    % Méthode d'Adams-Bashforth et Adams-Moulton
    for i = 3:n
        % Adams-Bashforth
        y_pred = y(:,i) + (h / 12) * (23 * feval(f, t(i), y(:,i)) ...
                - 16 * feval(f, t(i - 1), y(:,i-1)) + 5 * feval(f, t(i -2), y(:,i-2)));

        % Adams-Moulton
        y(:,i+1) = y(:,i) + (h / 24) * (9 * feval(f, t(i+1), y_pred) ...
                + 19 * feval(f, t(i), y(:,i)) - 5 * feval(f, t(i - 1), y(:,i-1)) ...
                + feval(f, t(i - 2), y(:,i-2)));
    end


end
