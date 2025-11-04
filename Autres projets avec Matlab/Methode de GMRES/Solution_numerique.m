% Il faut ce fichier avant le fichier(TP6_Q4_b.m)
clear all;
clc;
%% Question 4-a
m=200;

rng('default')
A=2*eye(m)+0.5*randn(m)/sqrt(m); 
b=ones(m,1);
tol=1e-20;

% la solution exacte
x_star=A\b;

    for k=1:24
            % [x, alpha_0, time ] = GMRES_Method(A, b, tol,iter)
        x_n=GMRES_Method(A, b, tol,k);
        Erreur(k)=norm(x_star-x_n,2);

    end

semilogy(Erreur,'r*','DisplayName','Erreur')
xlabel('n')
ylabel('Norme des erreurs')
title('Graphique de la convergence de la Methode GMRES')
grid on
legend('show')
%% Determination de C et D

n_b=1:24;

Erreur_Inf=norm(x_star-GMRES_Method(A, b, tol,2),2);
Erreur_Sup=norm(x_star-GMRES_Method(A, b, tol,24),2);
D=exp((log(Erreur_Inf)-log(Erreur_Sup))/(24-2));
C=exp(log(Erreur_Inf)+2*log(D));
    
        % Calcul de CD^(-n)

     CD_n=C*D.^(-n_b);
hold on
plot(n_b,CD_n,'b--','DisplayName','CD^{-n}')
hold off
















