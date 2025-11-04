function [Q_arnoldi, H_arnoldi] = arnoldi(A, Q, k)
    % k:nombre d'iteration
   
        V=A*Q(:,k);
    
        for i= 1:k     
            H_arnoldi(i) = Q(:,i)' * V;
            V = V - H_arnoldi(i) * Q(:,i);
        end
        
       H_arnoldi(k+1) = norm(V);
       Q_arnoldi= V / H_arnoldi(k+1);
    
end
