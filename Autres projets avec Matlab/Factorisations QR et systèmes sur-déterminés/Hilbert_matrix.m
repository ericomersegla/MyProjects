%TP4
clear all;
clc;

%Question 3

for n = 2:12
    A = Hilbert_matrix(n+10,n);
    x = ones(n,1);
    b = A*x;
    [q, r] = Trans_Householder(A);
    x_hat =r\(q'*b);
    delta = x-x_hat;
    res= b-A*x_hat;
    Err_rel(n) = norm(delta,2)/norm(x,2);
    Res_rel(n) = norm(res,2)/norm(b,2);
    kappa(n) = norm(A,2)*norm(pinv(A),2);
end
figure(1)
plot(2:12,Err_rel(2:12),'b-*')
hold on
plot(2:12,Res_rel(2:12),'r-o')
hold on
xlabel('n');
ylabel("erreurs et résidus relatufs");
title('Graphe des erreurs et résidus');
legend('Erreurs relatives', 'Résidus relatifs');
grid on

figure(2)
plot(2:12,kappa(2:12),'b-*')
xlabel('n');
ylabel("\kappa_2(A)");
title('Comportement du \kappa_2(A) en fonction de n');
grid on





function H = Hilbert_matrix(m, n)
    H = zeros(m, n);
    for i = 1:m
        for j = 1:n
            H(i, j) = 1 / (i + j - 1);
        end
    end
end



