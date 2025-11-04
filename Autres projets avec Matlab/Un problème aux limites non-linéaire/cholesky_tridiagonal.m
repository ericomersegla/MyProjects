% Vect_diag= [2, 2, 2];
% Vect_diag_sec= [-1, -1];
% [a,b]=cholesky(Vect_diag, Vect_diag_sec);

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
