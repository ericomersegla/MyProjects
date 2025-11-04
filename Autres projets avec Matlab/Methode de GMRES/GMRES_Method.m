
function [x, alpha_0, time ] = GMRES_Method(A, b, tol,iter)
    % iter ieteration maximale
    Q=zeros(length(b),iter);
    
    time=zeros(iter, 1);
    alpha_0=[];
    
    e_1=zeros(iter+1, 1);
    e_1(1)=1;
     

    c=zeros(iter, 1);
    s=zeros(iter, 1);
    
    Erreur=zeros(24, 1);
    % formation du vecteur de erreurs
  r = b;
  b_norm = norm(b,2);
  Erreur = norm(r,2) / b_norm;
    
 % formation de Q
  r_norm = norm(r);
  Q(:,1) = r / r_norm;
  
  alpha_1 = r_norm * e_1;
  
  for i = 1:iter
    tic; 

    % [Q_arnoldi, H_arnoldi] = arnoldi(A, Q, k)
    [Q(:, i+1),H(1:i+1, i) ] = arnoldi(A, Q, i);

    % [H_givens,c_given,s_given] = Rot_Givens(H,c,s,n)
    [H(1:i+1, i), c(i), s(i)] = Rot_Givens(H(1:i+1,i), c, s, i);

    time(i) = toc; % arreter le compteur de temps
    
    % Update the residual vector
    alpha_1(i + 1) = -s(i) * alpha_1(i);
    alpha_1(i)     = c(i) * alpha_1(i);
    Erreur       = abs(alpha_1(i + 1)) / b_norm;

    % Save the error
    alpha_0 = [alpha_0; Erreur];

    if (Erreur <= tol)
      break;
    end
  end
  
  % Calculate the result
  y = H(1:i, 1:i) \ alpha_1(1:i);
  x = Q(:, 1:i) * y;
end
   



