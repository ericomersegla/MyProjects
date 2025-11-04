

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
