%% personn
 clear all
 clc

a=0;
b=1;
n=100;
a2=-1*ones(n-2,1);
% n=(10^5 10^4 10^3 10^2 10^1 1);
y0=initialey0(a,b,n);
[a1,b1]=bcalcule(a,b,y0,n);
[delta1,z,y,x]=choleskysolve(a1,a2,b1);

%% fonction calcule de b et a1
function [a1,b1]=bcalcule(a,b,y0,n)
[x0,h]=xh(a,b,n);
a1=zeros(n-1,1);
b1=zeros(n-1,1);
a1(1)=2+h^2*exp(y0(1));
a1(n-1)=2+h^2*exp(y0(n-1));
b1(1)=-(2*y0(1)-y0(2)+h^2*(exp(y0(1))+2-exp(x0(1)^2)));
b1(n-1)=-(-y0(n-2)+2*y0(n-1)+h^2*(exp(y0(n-1))+2-exp(x0(n-1)^2)));
for j=2:n-2
b1(j)=-(-y0(j-1)+2*y0(j)-y0(j+1)+h^2*(exp(y0(j))+2-exp(x0(j)^2)));
a1(j)=2+h^2*exp(y0(j));
end
end
%% fonction initiale
function [y0]=initialey0(a,b,n)
y0=zeros(n-1,1);
[x0,~]=xh(a,b,n);
p=@(x) ((1/6)*(1-exp(1))*x^3+x^2+(1/6)*(-7+exp(1))*x);
for i=1:n-1
y0(i)=p(x0(i));
end
end
%% fonction calcule de x et h
function [x0,h]=xh(a,b,n)
h=(b-a)/n;
x0=zeros(n-1,1);
for i=1:n-1
x0(i)=a+i*h;
end
end

%% 2- fonction de cholesky modifiée
function [t,z,y,x]=choleskysolve(a1,a2,b)
n=size(a1,1);
x=zeros(n,1);
y=zeros(n-1,1);
z=zeros(n,1);
t=zeros(n,1);
x(1)=sqrt(a1(1));
z(1)=b(1)/x(1);
for i=2:n
y(i-1)=a2(i-1)/x(i-1);% diagonale secondaire en dessous de L
x(i)=sqrt(a1(i)-y(i-1)^2);% diagonale principale de L
if conj(x(i))~=x(i)%conj(x(i))~=x(i) %L(i,i) n'est pas réel
disp('matrice n`est pas HDP')
x=[];
return
end
if abs(x(i))<1e-16 %x(i) zéro
disp('matrice n`est pas HDP')
x=[];
return
end
z(i)=(b(i)-y(i-1)*z(i-1))/x(i);%substition avant
t(n)=z(n)/x(n);
for j=n-1:-1:1
t(j)=(z(j)-y(j)*t(j+1))/x(j);%substition arriere
end
end
end
