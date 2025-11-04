%TP4
clear all;
clc;

%Question 3

Y=2:12;

for i=2:8
    disp(i)
%Methode A Cholesky 
    A=Hilbert_matrix(i+10,i);
    B=A'*A;
    L=Fact_Cholesky(B); 
    Q=A*pinv(L'); 
    QT=Q';
    Q_2=QT*Q;
    z=size(Q_2);
    z=int64(z);
    I = eye( int64(z(2)), int64(z(2)) );  
    logerr_chol(i)=-log10(norm(I-QT*Q,2) );
end    



for i=2:12
    %Methode B Gram-Schmidt classique
    A=Hilbert_matrix(i+10,i);
    [Q,R]=Gram_Schmidt_Class(A); 

    QT=Q';
    Q_2=QT*Q;
    z=size(Q_2);
    z=int64(z);
    I = eye( int64(z(2)), int64( z(2) ) );
    logerr_gsc(i)=-log10(norm(I-QT*Q,2) ); 


%Methode C Gram-Schmidt modifiée
    [Q,R]=Gram_Schmidt_Mod(A); 
    QT=Q';
    logerr_mgs(i)=-log10(norm(I-QT*Q,2) );


  %Methode D Rotations de Givens  
    [Q,R]=Rot_Givens(A); 
    QT=Q';
    Q_2=QT*Q;
    z=size(Q_2);
    z=int64(z);
    I = eye( int64(z(2)), int64( z(2) ) );
    logerr_givens(i)=-log10(norm(I-QT*Q,2) );  

%Methode E Transformations de Householder
     x=ones(1,i);
     x=x';
     b=A*x;
     [Q,R]=Trans_Householder(A);
     QT=Q';
     Q_2=QT*Q;
     z=size(Q_2);
     z=int64(z);
     I = eye( int64(z(2)), int64( z(2) ) );
     logerr_hholder(i)=-log10(norm(I-QT*Q,2) );

end      
figure(1)

plot(2:8,logerr_chol(2:8),'r-*')
hold on;
plot(Y,logerr_gsc(2:12),'b-o')
hold on;
plot(Y,logerr_mgs(2:12),'k-o')
hold on;
plot(Y,logerr_givens(2:12),'r-o')
hold on;
plot(Y,logerr_hholder(2:12),'b-*')
legend('Cholesky','GSC' ,'GSM','Givens','HouseHolder' ),
hold on;
xlabel('n');
ylabel('-log_{10}(|| I_n - Q^T Q ||_2)');
title("Graphe des erreurs d'orthogonalisation en fonction des différentes méthodes");
grid on




function H = Hilbert_matrix(m, n)
    H = zeros(m, n);
    for i = 1:m
        for j = 1:n
            H(i, j) = 1 / (i + j - 1);
        end
    end
end
