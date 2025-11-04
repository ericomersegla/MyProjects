clear all; clc; 
%----------------------------Données pour l'option d'achat------------------------------------
s0=80;  %le prix de l'actif en début de période de l'entreprise spécialisée dans les voitures électriques
K_C=90; %le prix d'exercice
p_1=10;%la prime d'acquisition payée
kt=50; %nombre de simulation
N=Les_Instances(); %maturité de l'actif qui peut être les dates de retraites ou 1 mois avant les retraites à payer.
r=0.2; sigma_1=0.1; h_1=1.5; 
a=zeros(kt,N+1); %matrice de zero à utiliser dans la deuxième methode
Z=VA(kt,N); %v.a gaussienne

%----------------------------Données pour l'option de vente------------------------------------
x0=55;  %le prix de l'actif en début de période de l'entreprise spécialisée dans les hydrocarbures
K_P=100; %le prix d'exercice
p_2=11;%la prime d'acquisition payée
mu=0.3; sigma_2=0.5; h_2=0.75; 

%-------------------Première Methode: Methode du chainage arrière------------------------
St1=zeros(kt,N); %Inirialisation de la matrice des prix de l'option d'achat
Xt1=zeros(kt,N); %Inirialisation de la matrice des prix de l'option de vente
J1=zeros(kt,N+1); %Inirialisation de Matrices des coûts
JS1=zeros(kt,N+1); %Inirialisation de Matrices des coûts pour l'option d'achat
JX1=zeros(kt,N+1); %Inirialisation de Matrices des coûts pour l'option de vente

Jd1=zeros(kt,N); %Inirialisation de Matrices des coûts optimaux
Ta1=zeros(kt,N); %Initialisation de la matrice temps d'arrêt
St1(:,N)=s0*exp((r-sigma_1^2/2)*h_1+sigma_1*sqrt(h_1)*Z(:,N));
Xt1(:,N)=x0*exp((mu-sigma_2^2/2)*h_2+sigma_2*sqrt(h_2)*Z(:,N));
Yt1=St1+Xt1;
        %----------A  N, on a : ------------------------
J1(:,N)= max(0,St1(:,N)-K_C-p_1)+ max(0,K_P-p_2-Xt1(:,N)) ;  %Matrices des coûts
JS1(:,N)= max(0,St1(:,N)-K_C-p_1);
JX1(:,N)= max(0,K_P-p_2-Xt1(:,N)) ; 

Jd1(:,N)=max(0,St1(:,N)-K_C-p_1)+ max(0,K_P-p_2-Xt1(:,N)); %Matrices des coûts optimaux

for k=N-1:-1:1
    for j=1:kt
        St1(j,k)=s0*exp(((r-sigma_1^2)/2)*h_1+sigma_1*sqrt(h_1)*Z(j,k));
        Xt1(j,k)=x0*exp(((mu-sigma_2^2)/2)*h_2+sigma_2*sqrt(h_2)*Z(j,k));
        
        JS1(j,k)=max(St1(j,k)-K_C-p_1,exp(-r*h_1)*JS1(j,k+1));
        JX1(j,k)=max(K_P-p_2-Xt1(j,k),exp(-mu*h_2)*JX1(j,k+1));

        J1(j,k)=max(St1(j,k)-K_C-p_1,exp(-r*h_1)*JS1(j,k+1))+ max(K_P-p_2-Xt1(j,k),exp(-mu*h_2)*JX1(j,k+1))
    end
end
Jd1=J1;
D_1=[]; D_2=[];
for k=1:kt
    A1=max(Jd1(k,:));
    D_1(k)=A1;
    Index=find(Jd1(k,:)~=A1);
    Index1=find(Jd1(k,:)==A1);
    Jd1(k,Index)=0;
    D_2(k)=Index1;
    if A1==0
    else
        Ta1(k,Index1)=1; % Matrice des temps d'arret.
    end
end

% Calcul des indicateuers de position et de dispersion pour les coûts de chaque trajectoires
Moy_1=zeros(2,kt); ecart_type_1=zeros(2,kt); %initialisation de la matrice des couts
Moy_1(1,:)=1:kt; ecart_type_1(1,:)=1:kt;
for k=1:kt 
    Moy_1(2,k)=mean(J1(k,:)); %Calcul de la moyenne sur chaque trajectoire
    ecart_type_1(2,k)=std((J1(k,:))); %Calcul de l'ecart type sur chaque trajectoire
end

%-----------------------DeuxièmeMethode: Méthode Longstaff et Schwart--------------------
St2=trajectoire(s0,kt,N,Z,r,sigma_1,h_1) %matrice de prix du sous-jacent de l'option d'achat en fonction des trajectoire
Xt2=trajectoire(x0,kt,N,Z,mu,sigma_2,h_2) %matrice de prix du sous-jacent de l'optionde vente en fonction des trajectoire

JS2=zeros(kt,N+1); %Inirialisation de Matrices des coûts pour l'option d'achat
JX2=zeros(kt,N+1); %Inirialisation de Matrices des coûts pour l'option de vente
J2=zeros(kt,N+1); % Inirialisation de Matrices Matrices des coûts

Jd2=zeros(kt,N); %Matrices des coûts optimaux
Ta2=zeros(kt,N); %Initialisation de la matrice temps d'arrêt

%----------Au temps N on a : ------------------------
JS2(:,N+1)= max(a(:,N+1),St2(:,N+1)-K_C-p_1);
JX2(:,N+1)= max(a(:,N+1),K_P-p_2-Xt2(:,N+1)) ; 


J2(:,N+1)= JS2(:,N+1)+ JX2(:,N+1) ;  %Matrices des coûts à N
Jd2(:,N+1)=J2(:,N+1);  %Matrices des coût optimaux à N


        %----------Cas de l'option d'achat----------------

for k=N:-1:2
    t=0; X_C=[]; Y_C=[]; Tr_C=[];
    for j=1:kt
        if St2(j,k)<(K_C+p_1) % Cette condition détermine si l'option est dans la monnaie
            t=t+1;
            X_C(t)=St2(j,k);
            Y_C(t)=JS2(j,k+1)*exp(-r*h_1);
            JS2(t,k)=max(a(j,k), St2(j,k)-K_C-p_1);
        end
    end
    b_1=polyfit(X_C,Y_C,2); %Evaluation des coefficients de notre polynome
    Jc_1=polyval(b_1,X_C); %Cout espéré si on exerce l'option à t=k+1
    for j=1:kt
        if JS2(j,k)<Jc_1(1,j)
           JS2(j,k)=0 
        end
        [m n]=size(Jc_1);
        if j==n
               break
        end
    end
end

        %----------Cas de l'option de vente---------------

for k=N:-1:2
    t=0; X_P=[]; Y_P=[]; Tr_P=[];
    for j=1:kt
        if (Xt2(j,k)+p_2)<K_P % Cette condition détermine si l'option est dans la monnaie
            t=t+1;
            X_P(t)=Xt2(j,k);
            Y_P(t)=JX2(j,k+1)*exp(-mu*h_2);
            JX2(t,k)=max(a(j,k), K_P-p_2-Xt2(j,k));
        end
    end
    b_2=polyfit(X_P,Y_P,2); %Evaluation des coefficients de notre polynome
    Jc_2=polyval(b_2,X_P); %Cout espéré si on exerce l'option à t=k+1
    for j=1:kt
        if JX2(j,k)<Jc_2(1,j)
           JX2(j,k)=0 
        end
        [m n]=size(Jc_2);
        if j==n
               break
        end
    end
end

J2=JS2+JX2;

Jd2=J2; E_1=[]; E_2=[];


for k=1:kt
    A2=max(Jd2(k,:))
    E_1(k)=A2
    Index=find(Jd2(k,:)~=A2)
    Index1=find(Jd2(k,:)==A2)
    Jd2(k,Index)=0;
    if length(Index1)==1
        E_2(k)=Index1
    end
    if A2==0
    else
        Ta2(k,Index1)=1; % Matrice des temps d'arret
    end
end
J2(:,1)=[1:kt];
Ta2(:,1)=[1:kt];

% Calcul des indicateuers de position et de dispersion pour les coûts de chaque trajectoires
Moy_2=zeros(2,kt); ecart_type_2=zeros(2,kt); %initialisation de la matrice des couts
Moy_2(1,:)=1:kt; ecart_type_2(1,:)=1:kt;
for k=1:kt 
    Moy_2(2,k)=mean(J2(k,:)); %Calcul de la moyenne sur chaque trajectoire
    ecart_type_2(2,k)=std((J2(k,:))); %Calcul de l'ecart type sur chaque trajectoire
end

%1-  K_C+p_1<St2(j,k) & Xt2(j,k)<(K_P-p_2) ===== 0

%2-  K_C+p_1<St2(j,k) &(K_P-p_2)<Xt2(j,k) =======%max(St2(j,k)-K_C-p_1,Y_C(t))

%3-  St2(j,k)<K_C+p_1 & Xt2(j,k)<(K_P-p_2)===========%max(K_P-p_2-X_P(t),Y_P(t))
  
%4-  St2(j,k)< (K_C+p_1) & (K_P-p_2)<Xt2(j,k) ========= max(St2(j,k)-K_C-p_1,Y_C(t)) + max(K_P-p_2-X_P(t),Y_P(t))
%-------------pris en compte dans les algorithmes


%Analyses des résultats des deux méthodes
figure(1)
plot(D_1,'red','linewidth',2)
hold on
plot(E_1,'blue','linewidth',2)
grid on
xlabel('Trajectoires');
ylabel('Coûts Optimaux');
legend('Chaînage arrière','Longstaff et Schwartz')

figure(2)
plot(D_2,'red','linewidth',2)
hold on
plot(E_2,'blue','linewidth',2)
grid on
xlabel('Trajectoires');
ylabel('Temps d arrêt');
legend('Chaînage arrière','Longstaff et Schwartz')

%Calcul des erreurs sur la moyenne des coûts pour chaque méthodes
for k=1:kt
    Marge_erreur(k)=abs(Moy_1(2,k)-Moy_2(2,k))/Moy_1(2,k);
      
end
figure(3)
plot(Marge_erreur,'red','linewidth',2)
grid on
xlabel('Trajectoires');
ylabel('Erreurs Moyenne');
title('Erreur Relative')

% figure(4)
% plot(ecart_type_1,'red','Marker','*')
% grid on
% xlabel('Trajectoires');
% ylabel('Ecart type');
% legend('Chaînage arrière');

% figure(5)
% plot(ecart_type_2,'blue','Marker','*')
% grid on
% xlabel('Trajectoires');
% ylabel('Ecart type');
% legend('Longstaff et Schwartz');


% % verification de la trajectoire des moyenne "Moyenne de LS est
% % generalement en dessous dela moyenne de chainage arrière"
% figure(6)
% plot(Moy_1,'red','linewidth',2)
% hold on
% plot(Moy_2,'blue','linewidth',2)
% grid on
% xlabel('Trajectoires');
% ylabel('Moyenne');
% title('Moyenne des Coûts')
% legend('Chaînage arrière','Longstaff et Schwartz');

