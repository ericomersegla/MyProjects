%% 
clear all
clc


n=100;

a2=-1*ones(n-2,1);
a=0;
b=1;
% c1=2*ones(n,1);
% c2=-1*ones(n-1,1);
%[t,z,y,x]=choleskysolve(c1,c2,b);
[x0,h]=xh(a,b,n)
y0=initialey0(a,b,n)
[a1,b1]=bcalcule(a,b,y0,n);
[delta1,z,y,x]=choleskysolve(a1,a2,b1);
norm(delta1,Inf);
% for k=1:5
% n=10^k;
[y0,norm_residu,K,h]=newton(a,b,n);
%end

%% 3-b function de méthode de Newton

function [y0,norm_residu,K,h]=newton(a,b,n)
a2=-1*ones(n-2,1);
K=0;
y0=initialey0(a,b,n);
tol=1;
h=(b-a)/n;
norm_residu=0;
while tol>1e-10
[a1,b1]=bcalcule(a,b,y0,n);
[delta1,~,~,~]=choleskysolve(a1,a2,b1);

% Calcul de la norme infinie de la différence entre y(K) et y(K-1)
tol=norm(delta1, 'inf');
%la norme infinie du résidu final
norm_residu=norm(-b1,'inf');
y0=y0+delta1;
K=K+1;
end
y0=y0-delta1;
end
%% fonction calcule de b et a1
function [a1,b1]=bcalcule(a,b,y0,n)
[x0,h]=xh(a,b,n);
a1=zeros(n-1,1);
b1=zeros(n-1,1);
a1(1)=2+h^2*exp(y0(1));
a1(n-1)=2+h^2*exp(y0(n-1));
b1(1)=-(2*y0(1)-y0(2)+h^2*(exp(y0(1))+2-exp(x0(1))));
b1(n-1)=-(-y0(n-2)+2*y0(n-1)+h^2*(exp(y0(n-1))+2-exp(x0(n-1))));
for j=2:n-2
b1(j)=-(-y0(j-1)+2*y0(j)-y0(j+1)+h^2*(exp(y0(j))+2-exp(x0(j))));
a1(j)=2+h^2*exp(y0(j));
end
end
%% 3-a fonction initiale
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
