%% 
clear all
clc

a=0;
b=1;

for k=1:5
    n= 10^k;
    
Vect_diag_sec=-1*ones(n-2,1);

[x_ini,h]=Begin_x(a,b,n);
y0=Begin_y(a,b,n);
[alpha_1,beta_1]=factor_f(a,b,y0,n);
[Z_1,~,~,~]=cholesky_adapt(alpha_1,Vect_diag_sec,beta_1);

norm(Z_1,Inf);

[y0,K,h,norm_res,Iter_y]=newton(a,b,n)
end
for k=1:3 
    Err_iterr(k)=norm(Iter_y(:,5)-Iter_y(:,k),'inf'); 
    Err_iterr_1(k)=norm(Iter_y(:,5)-Iter_y(:,k+1),'inf'); 
 end

loglog(Err_iterr,Err_iterr_1);

p=polyfit(log(Err_iterr(1:2)), log(Err_iterr_1(1:2)),1);

disp(p);

%% 3-b function de méthode de Newton

function [y_ini,K,h,norm_res,Iter_y]=newton(a,b,n)
a2=-1*ones(n-2,1);
K=0;
y_ini=Begin_y(a,b,n);
tol=1;
h=(b-a)/n;
norm_res=0;
Iter_y=[y_ini]
while tol>1e-10
[a1,b1]=factor_f(a,b,y_ini,n);
[Z_1,~,~,~]=cholesky_adapt(a1,a2,b1);

% Calcul de la norme infinie de la différence entre y(K) et y(K-1)
tol=norm(Z_1, 'inf');
%la norme infinie du résidu final
norm_res=norm(-b1,'inf');
y_ini=y_ini+Z_1;
Iter_y=[Iter_y,y_ini];
K=K+1;
end
if K>=25
    fprintf(['Pas de convergence dans le nombre',...
         'd''itÈrations imparti\n ']);
    fprintf([' La valeur retourné a un rÈsidu',...
           relatif de e\n'],norm_res);
            fprintf(['La méthode a converge ‡ l''itération',...
            ' %i avec un rÈsidu %e\n'],K,norm_res);
end
end
%% fonction calcule de b et a1
function [alpha_1,beta_1]=factor_f(a,b,y_ini,n)
[x_ini,h]=Begin_x(a,b,n);
alpha_1=zeros(n-1,1);
beta_1=zeros(n-1,1);
alpha_1(1)=2+h^2*exp(y_ini(1));
alpha_1(n-1)=2+h^2*exp(y_ini(n-1));
beta_1(1)=-(2*y_ini(1)-y_ini(2)+h^2*(exp(y_ini(1))+2-exp(x_ini(1)^2)));
beta_1(n-1)=-(-y_ini(n-2)+2*y_ini(n-1)+h^2*(exp(y_ini(n-1))+2-exp(x_ini(n-1)^2)));
for j=2:n-2
beta_1(j)=-(-y_ini(j-1)+2*y_ini(j)-y_ini(j+1)+h^2*(exp(y_ini(j))+2-exp(x_ini(j)^2)));
alpha_1(j)=2+h^2*exp(y_ini(j));
end
end
%% 3-a fonction initiale
function [y_ini]=Begin_y(a,b,n)
y_ini=zeros(n-1,1);
[x_ini,~]=Begin_x(a,b,n);
fun=@(x) ((1/6)*(1-exp(1))*x^3+x^2+(1/6)*(-7+exp(1))*x);
for i=1:n-1
y_ini(i)=fun(x_ini(i));
end
end
%% fonction calcule de x et h
function [x_ini,h]=Begin_x(a,b,n)
h=(b-a)/n;
x_ini=zeros(n-1,1);
for i=1:n-1
x_ini(i)=a+i*h;
end
end

%% 2- fonction de cholesky modifiée
function [Z_1,Z_2,X_1,X_2]=cholesky_adapt(a1,a2,b)
n=size(a1,1);
X_2=zeros(n,1);
X_1=zeros(n-1,1);
Z_2=zeros(n,1);
Z_1=zeros(n,1);
X_2(1)=sqrt(a1(1));
Z_2(1)=b(1)/X_2(1);
for i=2:n
X_1(i-1)=a2(i-1)/X_2(i-1);% diagonale secondaire en dessous de L
X_2(i)=sqrt(a1(i)-X_1(i-1)^2);% diagonale principale de L
if conj(X_2(i))~=X_2(i)%conj(x(i))~=x(i) %L(i,i) n'est pas réel
disp('matrice n`est pas HDP')
X_2=[];
return
end
if abs(X_2(i))<1e-16 %x(i) zéro
disp('matrice n`est pas HDP')
X_2=[];
return
end
Z_2(i)=(b(i)-X_1(i-1)*Z_2(i-1))/X_2(i);%substition avant
Z_1(n)=Z_2(n)/X_2(n);
for j=n-1:-1:1
Z_1(j)=(Z_2(j)-X_1(j)*Z_1(j+1))/X_2(j);%substition arriere
end
end
end
