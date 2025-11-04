%% Region de stabilité
clear all
clc
%% Question 3)a)i)
method = 'AB';    
r = 3;            


switch method

  case 'euler'
     rho = @(z) -1 + z;
     sigma = @(z) 1;
     axisbox = [-3 1 -2 2];

  case 'AB'
     axisbox = [-3 1 -2 2];
     rho = @(z) (z-1) .* z.^(r-1);

     switch r
       case 1 
         sigma = @(z) ones(size(z));
       case 2 
         sigma = @(z) (3*z - 1)/2;
       case 3 
         sigma = @(z) ((23*z - 16) .* z +5) / 12;
       case 4 
         sigma = @(z) (((55*z - 59).*z +37).*z - 9) / 24;
       case 5 
         sigma = @(z) ((((1901*z - 2774).*z +2616).*z -1274).*z + 251) / 720;
       end

  case 'AM'
     % axisbox = [-7 1 -4 4];
     axisbox = [-7 1 -4 4];
     rho = @(z) (z-1) .* z.^(r-1);

     switch r
       case 1 
         sigma = @(z) (z + 1)/2;
       case 2 
         sigma = @(z) ((5*z + 8).*z -1) / 12;
       case 3 
         sigma = @(z) (((9*z + 19).*z -5).*z +1) / 24;
       case 4 
         sigma = @(z) ((((251*z + 646).*z - 264).*z + 106) -19).*z / 720;
       case 5 
         sigma = @(z) (((((475*z + 1427).*z - 798).*z + 482).*z ...
                       - 173).*z + 27)/ 1440;
       end
          
     axisbox = [-5 5 -5 5];

  otherwise
     disp(' *** methode absente')
  end
figure
plotBL(rho,sigma,axisbox)
xticks(-3:0.5:1);
xlabel('rho');
yticks(-2:0.5:2);
ylabel('sigma');
title('La région de stabilité absolue pour AB3');
grid on;
hold off;


function plotBL(rho,sigma,axisbox)
if nargin==3
    x_1 = axisbox(1); x_2 = axisbox(2); 
    y_1 = axisbox(3); y_2 = axisbox(4);
  else
    x_1 = -10; x_2 = 10; y_1 = -10; y_2 = 10;  % default values
end

clf
theta = linspace(0, 2*pi, 1001);
eitheta = exp(1i * theta);
z = rho(eitheta) ./ sigma(eitheta);
plot(real(z), imag(z))
box on
hold on

% plot axes:
plot([x_1 x_2],[0 0],'k')
plot([0 0],[y_1 y_2],'k')
axis([x_1 x_2 y_1 y_2])
set(gca,'FontSize',15)
hold off
end