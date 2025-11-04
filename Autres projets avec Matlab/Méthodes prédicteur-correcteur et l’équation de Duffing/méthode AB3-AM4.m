%% Définition des paramètres
clear all
clc

%% Définition des paramètres
g_m = 9.81;
l = 1;
gamma = 0.5;
omega = 1;
A = 0;

% Conditions initiales
alpha_0 = pi/16;
alpha_1 = 0;
n = 2.^(5:14);
h = 10./n;
t = 0:h:10;

y0 = [alpha_0; alpha_1]; % Correction des conditions initiales

% PARAMÈTRE DE LA SOLUTION
delta = (gamma^2) - 4 * (g_m/l);
r_1 = 0.5 * (-gamma + sqrt(delta));
r_2 = 0.5 * (-gamma - sqrt(delta));
c_1 = (r_2 * alpha_0 - alpha_1) / (r_2 - r_1);
c_2 = (alpha_1 - r_1 * alpha_0) / (r_2 - r_1);

 theta_1 = real(c_1 * exp(r_1 * t(end)) + c_2 * exp(r_2 * t(end)));
 theta_2 = real(r_1 * c_1 * exp(r_1 * t(end)) + r_2 * c_2 * exp(r_2 * t(end)));


% Initialisation de y pour stocker les résultats
y = zeros(2,length(n)); % Correction de la taille de y

for i = 1:length(n) % Correction de la boucle
    % Appel de la fonction Mafun avec les paramètres corrects
    [t_temp, y_temp] = Cauchy_AB3_AM4(@(t, y) Mafun(gamma, y, g_m, l), [0 10],y0, n(i));
    
    % Stocker la solution dans y
    y( :,i) = y_temp( :,end);
    

    % calcul des erreurs
    Error1(i)=sqrt(sum((y(1,i) - theta_1).^2));
    Error2(i)=sqrt(sum((y(2,i) - theta_2).^2));
    Error(i)=sqrt(sum((y(1,i) - theta_1).^2))+sqrt(sum((y(2,i) - theta_2).^2));
end

figure;
loglog(h, Error1, 'b-*', 'DisplayName', 'Erreur du predicteur AB3-AM4');
xlabel('h');
ylabel('Erreur');
title('Erreur avec la méthode du predicteur AB3-AM4 Tetha1');
legend('Location', 'Best');
grid on;

figure;
loglog(h, Error2, 'r-*', 'DisplayName', 'Erreur du predicteur AB3-AM4');
xlabel('h');
ylabel('Erreur');
title('Erreur avec la méthode du predicteur AB3-AM4 pour Tetha2');
legend('Location', 'Best');
grid on;


figure;
loglog(h, Error, 'b-*', 'DisplayName', 'Erreur du predicteur AB3-AM4');
xlabel('h');
ylabel('Erreur');
title('Somme des deux erreurs obtenues avec  predicteur AB3-AM4');
legend('Location', 'Best');
grid on;

hold off;

Ordre_con1=polyfit(log(h),log(Error1),1); % ordre de convergence
disp("L'ordre de convergence de la methode predicteur correcteur est de:" )
disp(Ordre_con1(1));
Ordre_con2=polyfit(log(h),log(Error2),1); % ordre de convergence
disp("L'ordre de convergence de la methode predicteur correcteur est de:" )
disp(Ordre_con2(1));
Ordre_con=polyfit(log(h),log(Error),1); % ordre de convergence
disp("L'ordre de convergence de la methode predicteur correcteur est de:" )
disp(Ordre_con(1));

function derive = Mafun(gamma, y, g, l) % Correction des paramètres
   derive = [y(2); -gamma * y(2) - (g * y(1))/l ]; % Ajout du terme A*sin(omega*t) 
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
