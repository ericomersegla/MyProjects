format shortG;
format compact;
t=0;
hmax=0.2;
tf=11.7;
gamma=[4*pi^2,4*pi^2,4*pi^2];
J=3;
N=20;
h=0.01;
x=[0,-(sqrt(3))/4,(sqrt(3))/4];
y=[1/2,-1/4,-1/4];
x1=[(3+sqrt(3)*cos(4*pi/9))/6,(-3+sqrt(3)*cos(4*pi/9))/6,2*sqrt(3)*cos(4*pi/9)/3];
y1=[(sqrt(3)*sin(4*pi/9))/6,(sqrt(3)*sin(4*pi/9))/6,(2*sqrt(3)*sin(4*pi/9))/3];
gamma1=[2,2,-1];
x_m= zeros(J,1000);
y_m= zeros(J, 1000);
x_m(:,1)=x1;
y_m(:,1)=y1;

[x_euler,y_euler]=euler(x,y,J,h,gamma);
[x_RK2,y_RK2]=RK2(x,y,J,h,gamma);
[x_RK4,y_RK4]=RK4(x,y,J,h,gamma);
[x_RK5,y_RK5]=RK5(x,y,J,h,gamma);
% [x_RKF,y_RKF]=RKF(x1,y1,J,hmax,tf,gamma1)
erreur_euler=erreur(x,y,N,J,gamma);
erreur_RK2=erreurRK2(x,y,N,J,gamma);
erreur_RK4=erreurRK4(x,y,N,J,gamma);
erreur_RK5=erreurRK5(x,y,N,J,gamma)
couleur=['r','b','g'];
%construction de la courbe de la méthode d'Euler
for i =1:J
plot(x_euler(i,:), y_euler(i,:),couleur(i),'LineWidth',2);
hold on;
end
title("Mouvement des tourbillons de la Méthode d'Euler");
legend('Tourbillon 1','Tourbillon 2','Tourbillon 3');
xlabel('x');
ylabel('y');
grid on;
hold off;
%construction de la courbe de la méthode de RK2
for i =1:J
plot(x_RK2(i,:), y_RK2(i,:), couleur(i),'LineWidth',1);
hold on;
end
title("Mouvement des tourbillons de la Méthode de RK2");
legend('Tourbillon 1','Tourbillon 2','Tourbillon 3');
xlabel('x');
ylabel('y');
grid on;
hold off;

%construction de la courbe de la méthode de RK4
for i=1:J
plot(x_RK4(i,:), y_RK4(i,:),couleur(i),'LineWidth',2);
hold on;
end
title("Mouvement des tourbillons de la Méthode de RK4");
legend('Tourbillon 1','Tourbillon 2','Tourbillon 3');
xlabel('x');
ylabel('y');
grid on;
hold off;

%construction de la courbe de la méthode de RK5
for i=1:J
plot(x_RK5(i,:), y_RK5(i,:),couleur(i),'LineWidth',2);
hold on;
end
title("Mouvement des tourbillons de la Méthode de RK5");
legend('Tourbillon 1','Tourbillon 2','Tourbillon 3');
xlabel('x');
ylabel('y');
grid on;
hold off;

%courbe d'erreur de la méthode d'Euler
H=1./2.^(1:20);
loglog(H, erreur_euler,'LineWidth',1)
hold on
title("erreur d'Euler")
xlabel('H');
ylabel('erreur_euler');
legend("L^2 Euler")
grid on;
% l'ordre de convergence de la méthode d'Euler
ordre=polyfit(H(2:7),erreur_euler(2:7),1)
hold off;
%courbe d'erreur de la méthode de RK2
H=1./2.^(1:20);
loglog(H, erreur_RK2,'LineWidth',1)
hold on;
title("erreur de RK2")
xlabel('H');
ylabel('erreur ');
legend("L^2 de RK2");
grid on;
hold off;
% l'ordre de convergence de la méthode de RK2
ordreRK2=polyfit(H(1:4),erreur_RK2(1:4),1)
%%%
%courbe d'erreur de la méthode de RK4
H=1./2.^(1:20);
loglog(H, erreur_RK4,"-o")
hold on;
title("erreur de RK4")
xlabel('H');
ylabel('erreur ');
%legend("L^2 de RK4");
grid on;
hold off;
% l'ordre de convergence de la méthode de RK4
ordreRK4=polyfit(H(3:6),erreur_RK4(3:6),1)

%courbe d'erreur de la méthode de RK5
H=1./2.^(1:20);
loglog(H, erreur_RK5,"o-")
hold on;
title("erreur de RK5")
xlabel('H');
ylabel('erreur ');
grid on;
hold off;
% l'ordre de convergence de la méthode de RK5
ordreRK5=polyfit(H(10:20),erreur_RK5(10:20),1)
%fionction d'erreur d'Euler

function erreur_euler=erreur(x,y,N,J,gamma)
x_erreur=zeros(J,1);
y_errur=zeros(J,1);
erreur_euler=zeros(1,N);
for i=1:N
h=1/2^i;
[x_euler,y_euler]=euler(x,y,J,h,gamma);
x_erreur(:,i)=x_euler(:,1)-x_euler(:,end);
y_erreur(:,i)=y_euler(:,1)-y_euler(:,end);
erreur_euler(i)=sum(norm([x_erreur(:,i),y_erreur(:,i)],2));
end
end
%Fonction d'erreur de RK2
function erreur_RK2=erreurRK2(x,y,N,J,gamma)
x_erreur=zeros(J,1);
y_erreur=zeros(J,1);
erreur_RK2=zeros(1,N);
for i=1:N
H=1/2^i;
[x_RK2,y_RK2]=RK2(x,y,J,H,gamma);
x_erreuR(:,i)=x_RK2(:,1)-x_RK2(:,end);
y_erreuR(:,i)=y_RK2(:,1)-y_RK2(:,end);
erreur_RK2(i)=sum(norm([x_erreuR(:,i),y_erreuR(:,i)],2));
end
end
%Fonction d'erreur de RK4
function erreur_RK4=erreurRK4(x,y,N,J,gamma)
x_erreuR4=zeros(J,1);
y_erreuR4=zeros(J,1);
erreur_RK4=zeros(1,N);
for i=1:N
H=1/2^i;
[x_RK4,y_RK4]=RK4(x,y,J,H,gamma);
x_erreuR4(:,i)=x_RK4(:,1)-x_RK4(:,end);
y_erreuR4(:,i)=y_RK4(:,1)-y_RK4(:,end);
erreur_RK4(i)=sum(norm([x_erreuR4(:,i),y_erreuR4(:,i)],2));
end
end
%Fonction d'erreur de RK5
function erreur_RK5=erreurRK5(x,y,N,J,gamma)
x_erreuR5=zeros(J,1);
y_erreuR5=zeros(J,1);
erreur_RK5=zeros(1,N);
for i=1:N
H=1/2^i;
[x_RK5,y_RK5]=RK5(x,y,J,H,gamma);
x_erreuR5(:,i)=x_RK5(:,1)-x_RK5(:,end);
y_erreuR5(:,i)=y_RK5(:,1)-y_RK5(:,end);
erreur_RK5(i)=sum(norm([x_erreuR5(:,i), y_erreuR5(:,i)],2));
end
end
%%Fonction de Runge kutta et Fortin
% function [x_RKF,y_RKF]=RKF(x1,y1,J,hmax,h,tf,k,gamma1)
% Atol=1e-13;
% Rtol=5e-10;
% while t<tf
% k1x=h*(fionctionx(x1,y1,gamma1,J));
% k1y=h*(fionctiony(x1,y1,gamma1,J));
% k2x=h*(fionctionx(0.25*k1x+x1,0.25*k1y+y1,gamma1,J));
% k2y=h*(fionctiony(0.25*k1x+x1,0.25*k1y+y1,gamma1,J));
% k3x=h*(fionctionx((3/32)*k1x+(9/32)*k2x+x1,(3/32)*k1y+(9/32)*k2y+y1,gamma1,J));
% k3y=h*(fionctiony((3/32)*k1x+(9/32)*k2x+x1,(3/32)*k1y+(9/32)*k2y+y1,gamma1,J));
% k4x=h*(fionctionx((1932/2197)*k1x-(7200/2197)*k2x+(7296/2197)*k3x+x1,(1932/2197)*k1y-(7200/2197)*k2y+(7296/2197)*k3y+y1,gamma1,J));
% k4y=h*(fionctiony((1932/2197)*k1x-(7200/2197)*k2x+(7296/2197)*k3x+x1,(1932/2197)*k1y-(7200/2197)*k2y+(7296/2197)*k3y+y1,gamma1,J));
% k5x=h*(fionctionx((439/216)*k1x-8*k2x+(3680/513)*k3x-(845/4104)*k4x+x1,(439/216)*k1y-8*k2y+(3680/513)*k3y-(845/4104)*k4y+y1,gamma1,J));
% k5y=h*(fionctiony((439/216)*k1x-8*k2x+(3680/513)*k3x-(845/4104)*k4x+x1,(439/216)*k1y-8*k2y+(3680/513)*k3y-(845/4104)*k4y+y1,gamma1,J));
% k6x=h*(fionctionx((-8/27)*k1x+2*k2x-(3544/2565)*k3x+(1859/4104)*k4x-(11/40)*k5x+x1,(-8/27)*k1y+2*k2y-(3544/2565)*k3y+(1859/4104)*k4y-(11/40)*k5y+y1,gamma1,J));
% k6y=h*(fionctiony((-8/27)*k1x+2*k2x-(3544/2565)*k3x+(1859/4104)*k4x-(11/40)*k5x+x1,(-8/27)*k1y+2*k2y-(3544/2565)*k3y+(1859/4104)*k4y-(11/40)*k5y+y1,gamma1,J));
% Ex=(1/360)*k1x-(128/4275)*k3x-(2197/75240)*k4x+(1/50)*k5x+(2/55)*k6x;
% Ey=(1/360)*k1y-(128/4275)*k3y-(2197/75240)*k4y+(1/50)*k5y+(2/55)*k6y;
% E=norm([Ex,Ey],"inf");
% tol=Atol+norm([x1,y1],"inf")*Rtol;
% x1=x1+(25/216)*k1x+(1408/2565)*k3x+(2197/4104)*k4x-(1/5)*k5x;
% y1=y1+(25/216)*k1y+(1408/2565)*k3y+(2197/4104)*k4y-(1/5)*k5y;
% beta=(tol/(2*E))^(1/5);
% h=min(beta*h,hmax);
% if E<=tol
% x_RKF(:,k)=x1;
% y_RKF(:,k)=y1;
% t=t+h;
% k=k+1;
% else
% h=h;
% RKF(x1,y1,J,hmax,h,tf,k,gamma1);
% 
% end
% end
% end
%%% Fonction de Runge kutta d'ordre 5
function [x_RK5,y_RK5]=RK5(x,y,J,h,gamma)
x_RK5 = zeros(J,101);
y_RK5= zeros(J, 101);
x_RK5(:,1)=x;
y_RK5(:,1)=y;
for k=1:100
k1x=h*(fionctionx(x_RK5(:,k),y_RK5(:,k),gamma,J));
k1y=h*(fionctiony(x_RK5(:,k),y_RK5(:,k),gamma,J));
k2x=h*(fionctionx(0.25*k1x+x_RK5(:,k),0.25*k1y+y_RK5(:,k),gamma,J));
k2y=h*(fionctiony(0.25*k1x+x_RK5(:,k),0.25*k1y+y_RK5(:,k),gamma,J));
k3x=h*(fionctionx((3/32)*k1x+(9/32)*k2x+x_RK5(:,k),(3/32)*k1y+(9/32)*k2y+y_RK5(:,k),gamma,J));
k3y=h*(fionctiony((3/32)*k1x+(9/32)*k2x+x_RK5(:,k),(3/32)*k1y+(9/32)*k2y+y_RK5(:,k),gamma,J));
k4x=h*(fionctionx((1932/2197)*k1x-(7200/2197)*k2x+(7296/2197)*k3x+x_RK5(:,k),(1932/2197)*k1y-(7200/2197)*k2y+(7296/2197)*k3y+y_RK5(:,k),gamma,J));
k4y=h*(fionctiony((1932/2197)*k1x-(7200/2197)*k2x+(7296/2197)*k3x+x_RK5(:,k),(1932/2197)*k1y-(7200/2197)*k2y+(7296/2197)*k3y+y_RK5(:,k),gamma,J));
k5x=h*(fionctionx((439/216)*k1x-8*k2x+(3680/513)*k3x-(845/4104)*k4x+x_RK5(:,k),(439/216)*k1y-8*k2y+(3680/513)*k3y-(845/4104)*k4y+y_RK5(:,k),gamma,J));
k5y=h*(fionctiony((439/216)*k1x-8*k2x+(3680/513)*k3x-(845/4104)*k4x+x_RK5(:,k),(439/216)*k1y-8*k2y+(3680/513)*k3y-(845/4104)*k4y+y_RK5(:,k),gamma,J));
k6x=h*(fionctionx((-8/27)*k1x+2*k2x-(3544/2565)*k3x+(1859/4104)*k4x-(11/40)*k5x+x_RK5(:,k),(-8/27)*k1y+2*k2y-(3544/2565)*k3y+(1859/4104)*k4y-(11/40)*k5y+y_RK5(:,k),gamma,J));
k6y=h*(fionctiony((-8/27)*k1x+2*k2x-(3544/2565)*k3x+(1859/4104)*k4x-(11/40)*k5x+x_RK5(:,k),(-8/27)*k1y+2*k2y-(3544/2565)*k3y+(1859/4104)*k4y-(11/40)*k5y+y_RK5(:,k),gamma,J));
x_RK5(:,k+1)=x_RK5(:,k)+(16/135)*k1x+(6656/12825)*k3x+(28561/56430)*k4x-(9/50)*k5x+(2/55)*k5x;
y_RK5(:,k+1)=y_RK5(:,k)+(16/135)*k1y+(6656/12825)*k3y+(28561/56430)*k4y-(9/50)*k5y+(2/55)*k5y;
end
end
%%% Fonction de Runge kutta d'ordre 4
function [x_RK4,y_RK4]=RK4(x,y,J,h,gamma)
x_RK4 = zeros(J,101);
y_RK4= zeros(J, 101);
x_RK4(:,1)=x;
y_RK4(:,1)=y;
for k=1:100
k1x=h*(fionctionx(x_RK4(:,k),y_RK4(:,k),gamma,J));
k1y=h*(fionctiony(x_RK4(:,k),y_RK4(:,k),gamma,J));
k2x=h*(fionctionx(0.25*k1x+x_RK4(:,k),0.25*k1y+y_RK4(:,k),gamma,J));
k2y=h*(fionctiony(0.25*k1x+x_RK4(:,k),0.25*k1y+y_RK4(:,k),gamma,J));
k3x=h*(fionctionx((3/32)*k1x+(9/32)*k2x+x_RK4(:,k),(3/32)*k1y+(9/32)*k2y+y_RK4(:,k),gamma,J));
k3y=h*(fionctiony((3/32)*k1x+(9/32)*k2x+x_RK4(:,k),(3/32)*k1y+(9/32)*k2y+y_RK4(:,k),gamma,J));
k4x=h*(fionctionx((1932/2197)*k1x-(7200/2197)*k2x+(7296/2197)*k3x+x_RK4(:,k),(1932/2197)*k1y-(7200/2197)*k2y+(7296/2197)*k3y+y_RK4(:,k),gamma,J));
k4y=h*(fionctiony((1932/2197)*k1x-(7200/2197)*k2x+(7296/2197)*k3x+x_RK4(:,k),(1932/2197)*k1y-(7200/2197)*k2y+(7296/2197)*k3y+y_RK4(:,k),gamma,J));
k5x=h*(fionctionx((439/216)*k1x-8*k2x+(3680/513)*k3x-(845/4104)*k4x+x_RK4(:,k),(439/216)*k1y-8*k2y+(3680/513)*k3y-(845/4104)*k4y+y_RK4(:,k),gamma,J));
k5y=h*(fionctiony((439/216)*k1x-8*k2x+(3680/513)*k3x-(845/4104)*k4x+x_RK4(:,k),(439/216)*k1y-8*k2y+(3680/513)*k3y-(845/4104)*k4y+y_RK4(:,k),gamma,J));
%k6x=h*(fionctionx((-8/27)*k1x+2*k2x-(3544/2565)*k3x+(1859/4104)*k4x-(11/40)*k5x+x_RK4(:,k),(-8/27)*k1y+2*k2y-(3544/2565)*k3y+(1859/4104)*k4y-(11/40)*k5y+y_RK4(:,k),gamma,J));
%k6y=h*(fionctiony((-8/27)*k1x+2*k2x-(3544/2565)*k3x+(1859/4104)*k4x-(11/40)*k5x+x_RK4(:,k),(-8/27)*k1y+2*k2y-(3544/2565)*k3y+(1859/4104)*k4y-(11/40)*k5y+y_RK4(:,k),gamma,J));
x_RK4(:,k+1)=x_RK4(:,k)+(25/216)*k1x+(1408/2565)*k3x+(2197/4104)*k4x-(1/5)*k5x;
y_RK4(:,k+1)=y_RK4(:,k)+(25/216)*k1y+(1408/2565)*k3y+(2197/4104)*k4y-(1/5)*k5y;
end
end

%%% Fonction de Runge kutha d'ordre 2
function [x_RK2,y_RK2]=RK2(x,y,J,h,gamma)
x_RK2 = zeros(J,101);
y_RK2= zeros(J, 101);
x_RK2(:,1)=x;
y_RK2(:,1)=y;
for k=1:100
x_RK2(:,k+1)=x_RK2(:,k)+0.5*h*(fionctionx(x_RK2(:,k),y_RK2(:,k),gamma,J)+fionctionx(x_RK2(:,k)+h*fionctionx(x_RK2(:,k),y_RK2(:,k),gamma,J), y_RK2(:,k)+h*fionctiony(x_RK2(:,k),y_RK2(:,k),gamma,J),gamma,J));
y_RK2(:,k+1)=y_RK2(:,k)+0.5*h*(fionctiony(x_RK2(:,k),y_RK2(:,k),gamma,J)+fionctiony(x_RK2(:,k)+h*fionctionx(x_RK2(:,k),y_RK2(:,k),gamma,J), y_RK2(:,k)+h*fionctiony(x_RK2(:,k),y_RK2(:,k),gamma,J),gamma,J));
end
end
%%%%%%%%
function [x_euler,y_euler]=euler(x,y,J,h,gamma)
x_euler = zeros(J,101);
y_euler = zeros(J, 101);
x_euler(:,1)=x;
y_euler(:,1)=y;
for k=1:100
x_euler(:,k+1)=x_euler(:,k)+h*fionctionx(x_euler(:,k),y_euler(:,k),gamma,J);
y_euler(:,k+1)=y_euler(:,k)+h*fionctiony(x_euler(:,k),y_euler(:,k),gamma,J);
end
end
%%%%%%%
function dx =fionctionx(x,y,gamma,J)
dx=zeros(J,1);
for i=1:J
for j=1:J
if j~=i
dx(i)=dx(i)-(gamma(j)*(y(i)-y(j))/((2*pi)*((x(i)-x(j))^2+(y(i)-y(j))^2)));
end
end
end
end
%%%%
function dy =fionctiony(x,y,gamma,J)
dy=zeros(J,1);
for i=1:J
for j=1:J
if j~=i
dy(i)=dy(i)+(gamma(j)*(x(i)-x(j))/((2*pi)*((x(i)-x(j))^2+(y(i)-y(j))^2)));
end
end
end
end
