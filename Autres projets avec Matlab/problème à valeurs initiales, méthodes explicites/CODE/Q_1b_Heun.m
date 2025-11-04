%% Refaire les mêmes calculs en utilisant la méthode de Heun (RK2).
x = [0, -sqrt(3)/4, sqrt(3)/4];
y = [1/2, -1/4, -1/4];

J = 3;
gamma = 4 * pi^2;
h = 0.01;
timesteps = 1 / h + 1;
z_val_heun = zeros(J + J, timesteps);
z_val_heun(:,1) = [transpose(x); transpose(y)]; % les premières valeurs sont les valeurs initiales

for i = 2:timesteps
    k1 = h * Res_Syst(z_val_heun(:, i - 1), J, gamma);
    k2 = h * Res_Syst(z_val_heun(:, i - 1) + k1, J, gamma);
    z_val_heun(:, i) = z_val_heun(:, i - 1) + 0.5 * (k1 + k2);
end

% Affichage des trajectoires
figure;
plot(z_val_heun(1,:), z_val_heun(4,:), "r-.");
hold on
plot(z_val_heun(2,:), z_val_heun(5,:), "b-.");
hold on
plot(z_val_heun(3,:), z_val_heun(6,:), "k-.");
xlabel("x")
ylabel("Y")
legend('Tourbillon 1', 'Tourbillon 2', 'Tourbillon 3')
title('Trajectoires la méthode de Heun ')
hold off

%% systeme dequation
function dzdt =Res_Syst(z,J, G)
    dzdt = zeros(J+J,1);

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
