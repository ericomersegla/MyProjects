function [Q,R] = Rot_Givens(A)
[m,n]=size(A);
Q = eye(m);
R = A;

for j = 1:n
    for i = m:-1:(j+1)
        G = eye(m);
        [c,s] = rotation(R(i-1,j),R(i,j));
        G(i-1,i-1)=c; G(i,i)=c; G(i-1,i)=s; G(i,i-1)=-s;
        Q = Q*G;
        R = G'*R;
    end
end
Q(:,1:m-1)=(-1)*Q(:,1:m-1);
R=(-1)*R;

function [c,s] = rotation(a,b)
if b == 0
   c = 1;
   s = 0;
 else
    if abs(b) > abs(a)
       r = -a / b;
       s = 1 / sqrt(1 + r^2);
       c = s*r;
    else
       r = -b / a;
       c = 1 / sqrt(1 + r^2);
       s = c*r;
    end
end