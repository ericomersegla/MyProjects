%% TP2 MAT6470_2024
clear all; 
clc;
%% Résolution par la méthode d'euler


% Positions initiales
x = [0, -sqrt(3)/4, sqrt(3)/4];
y = [1/2, -1/4, -1/4];

J = 3;
gamma = 4 * pi^2;
h = 0.01;
timesteps = 1 / h + 1;
z_val=zeros(J+J,timesteps);
z_val(:,1)= [transpose(x);transpose(y)] ;% les premieres valeurs sont les valeurs initiales, on a : 
for i= 2:timesteps
z_val(:,i)=z_val(:,i-1)+h*Res_Syst(z_val(:,i-1),J,gamma);
end

% Affichage des trajectoires
figure;
plot(z_val(1,:),z_val(4,:),"r-*");
hold on
plot(z_val(2,:),z_val(5,:),"b-*");
hold on
plot(z_val(3,:),z_val(6,:),"k-*");
xlabel("x")
ylabel("Y")
legend('Tourbillon 1','Tourbillon 2','Tourbillon 3')
title('Trajectoires des trois tourbillons par la méthode de Euler avec h=0.01 ')
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
