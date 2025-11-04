% test de la fonction
% f=@(t,y) y-t^2+1;
% tspan=[0, 10];
% n=2^5;
% y0=0.5;
%  [t, y_AM] =ffCauchy_AB3_AM4(f, tspan, y0, n);

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