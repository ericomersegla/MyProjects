function Xt=trajectoire(x0,kt,N,Z,mu,sigma,h)
Xt=zeros(kt,N+1); %initialisation du vecteurs prix de l'actif
Xt(:,1)=x0;
for k=1:N
    Xt(:,k+1)=x0*exp((mu-sigma^2/2)*(h)+sigma*sqrt(h)*Z(:,k));
end
    
