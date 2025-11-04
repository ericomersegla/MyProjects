%Question 2
t0=0;tf=11.7;
theta=(4*pi)/9;
x0=(3 + sqrt(3) * cos(theta)) / 6;y0=(sqrt(3) / 6) * sin(theta);
xx0=(-3 + sqrt(3) * cos(theta)) / 6;yy0=(sqrt(3) / 6) * sin(theta);
xxx0=(2 * sqrt(3) / 3) * cos(theta);yyy0=(2 * sqrt(3) / 3) * sin(theta);
h0=0.2;
Atol=1e-13;
Rtol=5e-10;

[t1,h1,u1,uu1,uuu1,v1,vv1,vvv1]=rk45(@fj1,@fj2,@fj3,@fj4,@fj5,@fj6,[0 12.6],x0,xx0,xxx0,y0,yy0,yyy0,10);
plot(u1,v1,'o');
hold on;
plot(uu1,vv1,'o');
hold on;
plot(uuu1,vvv1,'o');
xlabel('x(t)');ylabel('y(t)');
legend('RKF turbillon 1' , 'RKF turbillon 2','RKF turbillon 3');
title('RKF - trajectoires de trois tourbillions');
savefig('RKF trajectoires.fig')
close

plot(h1)
xlabel('t');ylabel('h(t)');
savefig('h.fig')
close

%Functions
function [t,u,uu,uuu,v,vv,vvv]=rk4(odef1,odef2,odef3,odef4,odef5,odef6,tspan,x0,xx0,xxx0,y0,yy0,yyy0,Nh)
tt=linspace(tspan(1),tspan(2),Nh+1);
h=(tspan(2)-tspan(1))/Nh;
u=x0;
uu=xx0;
uuu=xxx0;
v=y0;
vv=yy0;
vvv=yyy0;
for t=tt(1:end-1)
  x = u(end,:);
  xx = uu(end,:);
  xxx = uuu(end,:);
  y = v(end,:);
  yy = vv(end,:);
  yyy = vv(end,:);
  
  kx1=h*feval(odef1,t,x,xx,xxx,y,yy,yyy);
  kxx1=h*feval(odef2,t,x,xx,xxx,y,yy,yyy);
  kxxx1=h*feval(odef3,t,x,xx,xxx,y,yy,yyy);
  ky1=h*feval(odef4,t,x,xx,xxx,y,yy,yyy);
  kyy1=h*feval(odef5,t,x,xx,xxx,y,yy,yyy);
  kyyy1=h*feval(odef6,t,x,xx,xxx,y,yy,yyy);
  
  t1 = t + 0.25*h; 
  x1 = x + 0.25*kx1;
  xx1 = xx + 0.25*kxx1;
  xxx1 = xxx + 0.25*kxxx1;
  y1 = y + 0.25*ky1;
  yy1 = yy + 0.25*kyy1;
  yyy1 = yyy + 0.25*kyyy1;
  
  kx2=h*feval(odef1,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kxx2=h*feval(odef2,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kxxx2=h*feval(odef3,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  ky2=h*feval(odef4,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kyy2=h*feval(odef5,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kyyy2=h*feval(odef6,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  
  t2=t + (3/8)*h;
  x2 = x + (3/32)*kx1 + (9/32)*kx2;
  xx2 = xx + (3/32)*kxx1 + (9/32)*kxx2;
  xxx2 = xxx + (3/32)*kxxx1 + (9/32)*kxxx2;
  y2 = y + (3/32)*ky1 + (9/32)*ky2;
  yy2 = yy + (3/32)*kyy1 + (9/32)*kyy2;
  yyy2 = yyy + (3/32)*kyyy1 + (9/32)*kyyy2;
  
  kx3=h*feval(odef1,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kxx3=h*feval(odef2,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kxxx3=h*feval(odef3,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  ky3=h*feval(odef4,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kyy3=h*feval(odef5,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kyyy3=h*feval(odef6,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  
  t3=t + (12/13)*h; 
  x3 = x + (1932/2197)*kx1 - (7200/2197)*kx2 + (7296/2197)*kx3;
  xx3 = xx + (1932/2197)*kxx1 - (7200/2197)*kxx2 + (7296/2197)*kxx3;
  xxx3 = xxx + (1932/2197)*kxxx1 - (7200/2197)*kxxx2 + (7296/2197)*kxxx3;
  y3 = y + (1932/2197)*ky1 - (7200/2197)*ky2 + (7296/2197)*ky3;
  yy3 = yy + (1932/2197)*kyy1 - (7200/2197)*kyy2 + (7296/2197)*kyy3;
  yyy3 = yyy + (1932/2197)*kyyy1 - (7200/2197)*kyyy2 + (7296/2197)*kyyy3;
  
  kx4=h*feval(odef1,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kxx4=h*feval(odef2,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kxxx4=h*feval(odef3,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  ky4=h*feval(odef4,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kyy4=h*feval(odef5,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kyyy4=h*feval(odef6,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  
  t4=t+h;
  x4 = x + (439/216)*kx1 - 8*kx2 + (3680/513)*kx3 - (845/4104)*kx4;
  xx4 = xx + (439/216)*kxx1 - 8*kxx2 + (3680/513)*kxx3 - (845/4104)*kxx4;
  xxx4 = xxx + (439/216)*kxxx1 - 8*kxxx2 + (3680/513)*kxxx3 - (845/4104)*kxxx4;
  y4 = y + (439/216)*ky1 - 8*ky2 + (3680/513)*ky3 - (845/4104)*ky4;
  yy4 = yy + (439/216)*kyy1 - 8*kyy2 + (3680/513)*kyy3 - (845/4104)*kyy4;
  yyy4 = yyy + (439/216)*kyyy1 - 8*kyyy2 + (3680/513)*kyyy3 - (845/4104)*kyyy4;
  
  kx5=h*feval(odef1,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kxx5=h*feval(odef2,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kxxx5=h*feval(odef3,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  ky5=h*feval(odef4,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kyy5=h*feval(odef5,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kyyy5=h*feval(odef6,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  
  u = [u; u(end,:) + ((25/216)*kx1+(1408/2565)*kx3+(2197/4104)*kx4-(1/5)*kx5)];
  uu = [uu; uu(end,:) + ((25/216)*kxx1+(1408/2565)*kxx3+(2197/4104)*kxx4-(1/5)*kxx5)];
  uuu = [uuu; uuu(end,:) + ((25/216)*kxxx1+(1408/2565)*kxxx3+(2197/4104)*kxxx4-(1/5)*kxxx5)];
  v = [v; v(end,:) + ((25/216)*ky1+(1408/2565)*ky3+(2197/4104)*ky4-(1/5)*ky5)];
  vv = [vv; vv(end,:) + ((25/216)*kyy1+(1408/2565)*kyy3+(2197/4104)*kyy4-(1/5)*kyy5)];
  vvv = [vvv; vvv(end,:) + ((25/216)*kyyy1+(1408/2565)*kyyy3+(2197/4104)*kyyy4-(1/5)*kyyy5)];
end
end

function [t,u,uu,uuu,v,vv,vvv]=rk5(odef1,odef2,odef3,odef4,odef5,odef6,tspan,x0,xx0,xxx0,y0,yy0,yyy0,Nh)
tt=linspace(tspan(1),tspan(2),Nh+1);
h=(tspan(2)-tspan(1))/Nh;
u=x0;
uu=xx0;
uuu=xxx0;
v=y0;
vv=yy0;
vvv=yyy0;
for t=tt(1:end-1)
  x = u(end,:);
  xx = uu(end,:);
  xxx = uuu(end,:);
  y = v(end,:);
  yy = vv(end,:);
  yyy = vv(end,:);
  
  kx1=h*feval(odef1,t,x,xx,xxx,y,yy,yyy);
  kxx1=h*feval(odef2,t,x,xx,xxx,y,yy,yyy);
  kxxx1=h*feval(odef3,t,x,xx,xxx,y,yy,yyy);
  ky1=h*feval(odef4,t,x,xx,xxx,y,yy,yyy);
  kyy1=h*feval(odef5,t,x,xx,xxx,y,yy,yyy);
  kyyy1=h*feval(odef6,t,x,xx,xxx,y,yy,yyy);
  
  t1 = t + 0.25*h; 
  x1 = x + 0.25*kx1;
  xx1 = xx + 0.25*kxx1;
  xxx1 = xxx + 0.25*kxxx1;
  y1 = y + 0.25*ky1;
  yy1 = yy + 0.25*kyy1;
  yyy1 = yyy + 0.25*kyyy1;
  
  kx2=h*feval(odef1,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kxx2=h*feval(odef2,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kxxx2=h*feval(odef3,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  ky2=h*feval(odef4,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kyy2=h*feval(odef5,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kyyy2=h*feval(odef6,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  
  t2=t + (3/8)*h;
  x2 = x + (3/32)*kx1 + (9/32)*kx2;
  xx2 = xx + (3/32)*kxx1 + (9/32)*kxx2;
  xxx2 = xxx + (3/32)*kxxx1 + (9/32)*kxxx2;
  y2 = y + (3/32)*ky1 + (9/32)*ky2;
  yy2 = yy + (3/32)*kyy1 + (9/32)*kyy2;
  yyy2 = yyy + (3/32)*kyyy1 + (9/32)*kyyy2;
  
  kx3=h*feval(odef1,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kxx3=h*feval(odef2,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kxxx3=h*feval(odef3,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  ky3=h*feval(odef4,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kyy3=h*feval(odef5,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kyyy3=h*feval(odef6,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  
  t3=t + (12/13)*h; 
  x3 = x + (1932/2197)*kx1 - (7200/2197)*kx2 + (7296/2197)*kx3;
  xx3 = xx + (1932/2197)*kxx1 - (7200/2197)*kxx2 + (7296/2197)*kxx3;
  xxx3 = xxx + (1932/2197)*kxxx1 - (7200/2197)*kxxx2 + (7296/2197)*kxxx3;
  y3 = y + (1932/2197)*ky1 - (7200/2197)*ky2 + (7296/2197)*ky3;
  yy3 = yy + (1932/2197)*kyy1 - (7200/2197)*kyy2 + (7296/2197)*kyy3;
  yyy3 = yyy + (1932/2197)*kyyy1 - (7200/2197)*kyyy2 + (7296/2197)*kyyy3;
  
  kx4=h*feval(odef1,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kxx4=h*feval(odef2,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kxxx4=h*feval(odef3,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  ky4=h*feval(odef4,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kyy4=h*feval(odef5,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kyyy4=h*feval(odef6,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  
  t4=t+h;
  x4 = x + (439/216)*kx1 - 8*kx2 + (3680/513)*kx3 - (845/4104)*kx4;
  xx4 = xx + (439/216)*kxx1 - 8*kxx2 + (3680/513)*kxx3 - (845/4104)*kxx4;
  xxx4 = xxx + (439/216)*kxxx1 - 8*kxxx2 + (3680/513)*kxxx3 - (845/4104)*kxxx4;
  y4 = y + (439/216)*ky1 - 8*ky2 + (3680/513)*ky3 - (845/4104)*ky4;
  yy4 = yy + (439/216)*kyy1 - 8*kyy2 + (3680/513)*kyy3 - (845/4104)*kyy4;
  yyy4 = yyy + (439/216)*kyyy1 - 8*kyyy2 + (3680/513)*kyyy3 - (845/4104)*kyyy4;
  
  kx5=h*feval(odef1,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kxx5=h*feval(odef2,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kxxx5=h*feval(odef3,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  ky5=h*feval(odef4,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kyy5=h*feval(odef5,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kyyy5=h*feval(odef6,t4,x4,xx4,xxx4,y4,yy4,yyy4);
 
  t5=t+0.5*h;
  x5 = x - (8/27)*kx1 + 2*kx2 - (3544/2565)*kx3 + (1859/4104)*kx4 - (11/40)*kx5;
  xx5 = xx - (8/27)*kxx1 + 2*kxx2 - (3544/2565)*kxx3 + (1859/4104)*kxx4 - (11/40)*kxx5;
  xxx5 = xxx - (8/27)*kxxx1 + 2*kxxx2 - (3544/2565)*kxxx3 + (1859/4104)*kxxx4 - (11/40)*kxxx5;
  y5 = y - (8/27)*ky1 + 2*ky2 - (3544/2565)*ky3 + (1859/4104)*ky4 - (11/40)*ky5;
  yy5 = yy - (8/27)*kyy1 + 2*kyy2 - (3544/2565)*kyy3 + (1859/4104)*kyy4 - (11/40)*kyy5;
  yyy5 = yyy - (8/27)*kyyy1 + 2*kyyy2 - (3544/2565)*kyyy3 + (1859/4104)*kyyy4 - (11/40)*kyyy5;
  
  kx6=h*feval(odef1,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  kxx6=h*feval(odef2,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  kxxx6=h*feval(odef3,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  ky6=h*feval(odef4,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  kyy6=h*feval(odef5,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  kyyy6=h*feval(odef6,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  
  u = [u; u(end,:) + ((16/135)*kx1+(6656/12825)*kx3+(28561/56430)*kx4-(9/50)*kx5+(2/55)*kx6)];
  uu = [uu; uu(end,:) + ((16/135)*kxx1+(6656/12825)*kxx3+(28561/56430)*kxx4-(9/50)*kxx5+(2/55)*kxx6)];
  uuu = [uuu; uuu(end,:) + ((16/135)*kxxx1+(6656/12825)*kxxx3+(28561/56430)*kxxx4-(9/50)*kxxx5+(2/55)*kxxx6)];
  v = [v; v(end,:) + ((16/135)*ky1+(6656/12825)*ky3+(28561/56430)*ky4-(9/50)*ky5+(2/55)*ky6)];
  vv = [vv; vv(end,:) + ((16/135)*kyy1+(6656/12825)*kyy3+(28561/56430)*kyy4-(9/50)*kyy5+(2/55)*kyy6)];
  vvv = [vvv; vvv(end,:) + ((16/135)*kyyy1+(6656/12825)*kyyy3+(28561/56430)*kyyy4-(9/50)*kyyy5+(2/55)*kyyy6)];
end
end

function [t,h,u,uu,uuu,v,vv,vvv]=rk45(odef1,odef2,odef3,odef4,odef5,odef6,tspan,x0,xx0,xxx0,y0,yy0,yyy0,Nh)
h=0.2;
hmax=0.2;
Atol=1e-13;
Rtol=5e-10;
tt=linspace(tspan(1),tspan(2),Nh+1);
u=x0;
uu=xx0;
uuu=xxx0;
v=y0;
vv=yy0;
vvv=yyy0;
for t=tt(1:end-1)
  x = u(end,:);
  xx = uu(end,:);
  xxx = uuu(end,:);
  y = v(end,:);
  yy = vv(end,:);
  yyy = vv(end,:);
  
  kx1=h*feval(odef1,t,x,xx,xxx,y,yy,yyy);
  kxx1=h*feval(odef2,t,x,xx,xxx,y,yy,yyy);
  kxxx1=h*feval(odef3,t,x,xx,xxx,y,yy,yyy);
  ky1=h*feval(odef4,t,x,xx,xxx,y,yy,yyy);
  kyy1=h*feval(odef5,t,x,xx,xxx,y,yy,yyy);
  kyyy1=h*feval(odef6,t,x,xx,xxx,y,yy,yyy);
  
  t1 = t + 0.25*h; 
  x1 = x + 0.25*kx1;
  xx1 = xx + 0.25*kxx1;
  xxx1 = xxx + 0.25*kxxx1;
  y1 = y + 0.25*ky1;
  yy1 = yy + 0.25*kyy1;
  yyy1 = yyy + 0.25*kyyy1;
  
  kx2=h*feval(odef1,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kxx2=h*feval(odef2,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kxxx2=h*feval(odef3,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  ky2=h*feval(odef4,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kyy2=h*feval(odef5,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  kyyy2=h*feval(odef6,t1,x1,xx1,xxx1,y1,yy1,yyy1);
  
  t2=t + (3/8)*h;
  x2 = x + (3/32)*kx1 + (9/32)*kx2;
  xx2 = xx + (3/32)*kxx1 + (9/32)*kxx2;
  xxx2 = xxx + (3/32)*kxxx1 + (9/32)*kxxx2;
  y2 = y + (3/32)*ky1 + (9/32)*ky2;
  yy2 = yy + (3/32)*kyy1 + (9/32)*kyy2;
  yyy2 = yyy + (3/32)*kyyy1 + (9/32)*kyyy2;
  
  kx3=h*feval(odef1,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kxx3=h*feval(odef2,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kxxx3=h*feval(odef3,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  ky3=h*feval(odef4,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kyy3=h*feval(odef5,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  kyyy3=h*feval(odef6,t2,x2,xx2,xxx2,y2,yy2,yyy2);
  
  t3=t + (12/13)*h; 
  x3 = x + (1932/2197)*kx1 - (7200/2197)*kx2 + (7296/2197)*kx3;
  xx3 = xx + (1932/2197)*kxx1 - (7200/2197)*kxx2 + (7296/2197)*kxx3;
  xxx3 = xxx + (1932/2197)*kxxx1 - (7200/2197)*kxxx2 + (7296/2197)*kxxx3;
  y3 = y + (1932/2197)*ky1 - (7200/2197)*ky2 + (7296/2197)*ky3;
  yy3 = yy + (1932/2197)*kyy1 - (7200/2197)*kyy2 + (7296/2197)*kyy3;
  yyy3 = yyy + (1932/2197)*kyyy1 - (7200/2197)*kyyy2 + (7296/2197)*kyyy3;
  
  kx4=h*feval(odef1,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kxx4=h*feval(odef2,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kxxx4=h*feval(odef3,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  ky4=h*feval(odef4,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kyy4=h*feval(odef5,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  kyyy4=h*feval(odef6,t3,x3,xx3,xxx3,y3,yy3,yyy3);
  
  t4=t+h;
  x4 = x + (439/216)*kx1 - 8*kx2 + (3680/513)*kx3 - (845/4104)*kx4;
  xx4 = xx + (439/216)*kxx1 - 8*kxx2 + (3680/513)*kxx3 - (845/4104)*kxx4;
  xxx4 = xxx + (439/216)*kxxx1 - 8*kxxx2 + (3680/513)*kxxx3 - (845/4104)*kxxx4;
  y4 = y + (439/216)*ky1 - 8*ky2 + (3680/513)*ky3 - (845/4104)*ky4;
  yy4 = yy + (439/216)*kyy1 - 8*kyy2 + (3680/513)*kyy3 - (845/4104)*kyy4;
  yyy4 = yyy + (439/216)*kyyy1 - 8*kyyy2 + (3680/513)*kyyy3 - (845/4104)*kyyy4;
  
  kx5=h*feval(odef1,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kxx5=h*feval(odef2,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kxxx5=h*feval(odef3,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  ky5=h*feval(odef4,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kyy5=h*feval(odef5,t4,x4,xx4,xxx4,y4,yy4,yyy4);
  kyyy5=h*feval(odef6,t4,x4,xx4,xxx4,y4,yy4,yyy4);
 
  t5=t+0.5*h;
  x5 = x - (8/27)*kx1 + 2*kx2 - (3544/2565)*kx3 + (1859/4104)*kx4 - (11/40)*kx5;
  xx5 = xx - (8/27)*kxx1 + 2*kxx2 - (3544/2565)*kxx3 + (1859/4104)*kxx4 - (11/40)*kxx5;
  xxx5 = xxx - (8/27)*kxxx1 + 2*kxxx2 - (3544/2565)*kxxx3 + (1859/4104)*kxxx4 - (11/40)*kxxx5;
  y5 = y - (8/27)*ky1 + 2*ky2 - (3544/2565)*ky3 + (1859/4104)*ky4 - (11/40)*ky5;
  yy5 = yy - (8/27)*kyy1 + 2*kyy2 - (3544/2565)*kyy3 + (1859/4104)*kyy4 - (11/40)*kyy5;
  yyy5 = yyy - (8/27)*kyyy1 + 2*kyyy2 - (3544/2565)*kyyy3 + (1859/4104)*kyyy4 - (11/40)*kyyy5;
  
  kx6=h*feval(odef1,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  kxx6=h*feval(odef2,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  kxxx6=h*feval(odef3,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  ky6=h*feval(odef4,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  kyy6=h*feval(odef5,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  kyyy6=h*feval(odef6,t5,x5,xx5,xxx5,y5,yy5,yyy5);
  
  Ex=max(abs(((1/360)*kx1-(128/4275)*kx3-(2197/75240)*kx4+(1/50)*kx5+(2/55)*kx6)));
  Exx=max(abs(((1/360)*kxx1-(128/4275)*kxx3-(2197/75240)*kxx4+(1/50)*kxx5+(2/55)*kxx6)));
  Exxx=max(abs(((1/360)*kxxx1-(128/4275)*kxxx3-(2197/75240)*kxxx4+(1/50)*kxxx5+(2/55)*kxxx6)));
  Ey=max(abs(((1/360)*ky1-(128/4275)*ky3-(2197/75240)*ky4+(1/50)*ky5+(2/55)*ky6)));
  Eyy=max(abs(((1/360)*kyy1-(128/4275)*kyy3-(2197/75240)*kyy4+(1/50)*kyy5+(2/55)*kyy6)));
  Eyyy=max(abs(((1/360)*kyyy1-(128/4275)*kyyy3-(2197/75240)*kyyy4+(1/50)*kyyy5+(2/55)*kyyy6)));
  
  normErr=norm([Ex Ey;Exx Eyy;Exxx Eyyy]);
  tol=Atol+norm([x y;xx yy;xxx yyy])*Rtol;
  beta=(tol/normErr)^(1/5);
  h=h*beta
  
  u = [u; u(end,:) + ((16/135)*kx1+(6656/12825)*kx3+(28561/56430)*kx4-(9/50)*kx5+(2/55)*kx6)];
  uu = [uu; uu(end,:) + ((16/135)*kxx1+(6656/12825)*kxx3+(28561/56430)*kxx4-(9/50)*kxx5+(2/55)*kxx6)];
  uuu = [uuu; uuu(end,:) + ((16/135)*kxxx1+(6656/12825)*kxxx3+(28561/56430)*kxxx4-(9/50)*kxxx5+(2/55)*kxxx6)];
  v = [v; v(end,:) + ((16/135)*ky1+(6656/12825)*ky3+(28561/56430)*ky4-(9/50)*ky5+(2/55)*ky6)];
  vv = [vv; vv(end,:) + ((16/135)*kyy1+(6656/12825)*kyy3+(28561/56430)*kyy4-(9/50)*kyy5+(2/55)*kyy6)];
  vvv = [vvv; vvv(end,:) + ((16/135)*kyyy1+(6656/12825)*kyyy3+(28561/56430)*kyyy4-(9/50)*kyyy5+(2/55)*kyyy6)];
end
end

function fjx=fj1(t,x,xx,xxx,y,yy,yyy)
gamma1 = 2;
gamma2 = 2;
gamma3=-1;
fjx = -(gamma2/(2*pi))*(y-yy)/((x-xx)^2 + (y-yy)^2) - (gamma3/(2*pi))*(y-yyy)/((x-xxx)^2 + (y-yyy)^2); 
end

function fjxx=fj2(t,x,xx,xxx,y,yy,yyy)
gamma1 = 2;
gamma2 = 2;
gamma3=-1;
fjxx = -(gamma1/(2*pi))*(yy-y)/((xx-x)^2 + (yy-y)^2) - (gamma3/(2*pi))*(yy-yyy)/((xx-xxx)^2 + (yy-yyy)^2); 
end

function fjxxx=fj3(t,x,xx,xxx,y,yy,yyy)
gamma1 = 2;
gamma2 = 2;
gamma3=-1;
fjxxx = -(gamma1/(2*pi))*(yyy-y)/((xxx-x)^2 + (yyy-y)^2) - (gamma2/(2*pi))*(yyy-yy)/((xxx-xx)^2 + (yyy-yy)^2); 
end

function fjy=fj4(t,x,xx,xxx,y,yy,yyy)
gamma1 = 2;
gamma2 = 2;
gamma3=-1;
fjy = (gamma2/(2*pi))*(x-xx)/((x-xx)^2 + (y-yy)^2) + (gamma3/(2*pi))*(x-xxx)/((x-xxx)^2 + (y-yyy)^2);
end

function fjyy=fj5(t,x,xx,xxx,y,yy,yyy)
gamma1 = 2;
gamma2 = 2;
gamma3=-1;
fjyy=(gamma1/(2*pi))*(xx-x)/((xx-x)^2 + (yy-y)^2) + (gamma3/(2*pi))*(xx-xxx)/((xx-xxx)^2 + (yy-yyy)^2);
end

function fjyyy=fj6(t,x,xx,xxx,y,yy,yyy)
gamma1 =2;
gamma2 = 2;
gamma3=-1;
fjyyy=(gamma1/(2*pi))*(xxx-x)/((xxx-x)^2 + (yyy-y)^2) + (gamma2/(2*pi))*(xxx-xx)/((xxx-xx)^2 + (yyy-yy)^2); 
end