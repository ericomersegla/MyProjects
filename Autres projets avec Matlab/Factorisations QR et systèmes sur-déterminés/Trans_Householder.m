function [Q,R]=Trans_Householder(A)
    %Methode de Housesolder
    [m,n]=size(A); 
    Q=eye(m); 
    R=zeros(m,n);
    Al=A;
    for k=1:n
        e=zeros(1,m); e(1,k)=1;
        x=zeros(m,1);
        x(k:m,1)=Al(k:m,k);
        v=sign(x(k,1))*norm(x)*e'+x;
        v=v/norm(v);
        Q=Q-2*Q*v*v';
        Al=Al-2*v*v'*Al;
    end
    R=Al;
end


