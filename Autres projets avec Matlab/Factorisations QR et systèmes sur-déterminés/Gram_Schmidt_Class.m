function [Q,R]=Gram_Schmidt_Class(A)
[m,n]=size(A); Q=zeros(m,n); R=zeros(n,n); I=eye(m);

for j=1:n
    if j==1
       Q(:,j)=A(:,j);
    else
        Pp=I;
        for i=1:j-1
            R(i,j)=Q(:,i)'*A(:,j);
            P=Q(:,i)*Q(:,i)';
            Pp=Pp-P;
        end
        Q(:,j)=Pp*A(:,j);
    end
    R(j,j)=norm(Q(:,j));
    Q(:,j)=Q(:,j)/R(j,j);
end

