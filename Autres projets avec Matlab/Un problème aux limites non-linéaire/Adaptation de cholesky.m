%% 3-i Utilisons cholesky developpé en 2 pour resoudre le système linéaire à chaque itératif:
 clear all
 clc

% a=0;
% b=1;
% 
% n=100;
% h=1./n;
% 
% 
% Vect_diag_sec=-1*ones(n-2,1);
% [x_ini,y_ini]=Begin(a,b,n);
% [Vect_diag,x]=Factor(a,b,n);
% [Z_1,Z_2] =Resol_Sys(Vect_diag, Vect_diag_sec,x) ;


a = 0;
b = 1;
h_values = 10.^(-5:-1);
tol = 1;

for h = h_values
    n = round(1/h); % Assurer que n est un entier

    Vect_diag_sec = -1 * ones(n - 2, 1); 

    [x_ini, y_ini] = Begin(a, b, n);
    [Vect_diag, x] = Factor(a, b, n);
    [Z_1,Z_2] = Resol_Sys(Vect_diag, Vect_diag_sec, x);

    
    % disp(['Pour les valeurs de x = ', num2str(Z_1), ', Pour les valeurs de y = = ', num2str(Z_2)]);
end

%% Determination des paramétres

function [alpha,beta] =Factor(a, b,n) 
alpha=zeros(n-1,1);
beta=zeros(n-1,1);
h=(b-a)/n;
[x_ini,y_ini]=Begin(a,b,n);

alpha(1)=2+h^2*exp(y_ini(1));
alpha(n-1)=2+h^2*exp(y_ini(n-1));

beta(1)=-(2*y_ini(1)-y_ini(2)+h^2*(exp(y_ini(1))+2-exp(x_ini(1))));
beta(n-1)=-(-y_ini(n-2)+2*y_ini(n-1)+h^2*(exp(y_ini(n-1))+2-exp(x_ini(n-1))));

for k=2:n-2

beta(k)=-(-y_ini(k-1)+2*y_ini(k)-y_ini(k+1)+h^2*(exp(y_ini(k))+2-exp(x_ini(k))));
alpha(k)=2+h^2*exp(y_ini(k));

end

end

%% Initialisation
function [x_ini,y_ini]=Begin(a,b,n)
% a et b sont les bornes des intervalles et n le nombre de pas
y_ini=zeros(n-1,1);
x_ini=zeros(n-1,1);
h=(b-a)/n;
fun=@(y)((1/6)*(1-exp(1))*y^3+y^2+(1/6)*(-7+exp(1))*y);
for k=1:n-1
    x_ini(k)=a+k*h;
    y_ini(k)=fun(x_ini(k));
end

end

%% Resolution du système
function [Z_1,Z_2] = Resol_Sys(Vect_diag, Vect_diag_sec, x)
    n = length(Vect_diag);
    Z_1 =zeros(n,1);
    Z_2 =zeros(n,1);
    [X_1, X_2] = cholesky_tridiagonal(Vect_diag, Vect_diag_sec);
    Z_2(1) = x(1) / X_1(1);
    for k = 2:n
        Z_2(k) = (x(k) - X_2(k-1) * Z_2(k-1)) / X_1(k);
    end
    Z_1(n) = Z_2(n) / X_1(n);
    for i = n-1:-1:1
        Z_1(i) = (Z_2(i) - X_2(i) * Z_1(i+1)) / X_1(i);
    end
end


%% Adaptation de cholesky
function [X_1,X_2] = cholesky_tridiagonal(Vect_diag, Vect_diag_sec) %cholesky_tridiagonal
% Vect_diag: vecteur des éléments diagonaux de la matrice tridiagonale
% Vect_diag_sec: vecteur des éléments hors-diagonaux de la matrice tridiagonale
n=length(Vect_diag);
X_2 =zeros(n-1,1);
X_1 =zeros(n,1);
X_1(1)=sqrt(Vect_diag(1));

% Étape de décomposition de Cholesky
for k =2:n
    X_2(k-1)=Vect_diag_sec(k-1)/X_1(k-1);
    X_1(k)=sqrt(Vect_diag(k)-(X_2(k-1))^2); 
        
    % Mise à jour des éléments du vecteur du côté droit
    if conj(X_1(k))~= X_1(k)
        disp("Matrice n'est pas HDP")
        X_1=[];
        X_2=[];
           return
    end
    
    if abs(X_1(k))<1e-16
       disp("Matrice n'est pas HDP")
        X_1=[];
        X_2=[];
           return
    end
end
end
