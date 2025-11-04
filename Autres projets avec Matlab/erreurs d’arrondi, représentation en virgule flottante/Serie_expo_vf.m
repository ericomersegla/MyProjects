clear all;
clc
%% 
%2- Fonction pour calculer S(x;N)
% Valeurs de x
x_values = [-5, 5];

% Valeurs de N
N_values = 0:50;

% Calcul et stockage des erreurs relatives et de la borne superieure

% Initialisation des matrices pour stocker les résultats
Erreur_relative_1 = zeros(1, length(N_values));
Erreur_relative_2 = zeros(1, length(N_values));
Borne_Sup = zeros(1, length(N_values));
Borne_Sup2 = zeros(1, length(N_values));
S_value1=zeros(1, length(N_values));
S_value2=zeros(1, length(N_values));
 

%% 
% Tracé du graphique semilogy
figure;
semilogy(0:50, Erreur_relative(-5)); 
hold on;
semilogy(0:50, Erreur_relative(5));
hold on
semilogy(0:50, bornes(-5));
hold on
semilogy(0:50, bornes(5));
xlabel('N');
ylabel('Erreur relative');
legend('x = -5', 'x = 5', 'Borne Superieure théorique pour x = -5', 'Borne Superieure théorique pour x = 5');
title('Erreur relative de S(x;N) par rapport à exp(x)');



%%  3-Fonction pour calculer S(x;N) en précision simple
% Calcul des sommes partielles en ordre croissant et décroissant

%Pour x=-5
x_val = -5;
Sl1=single(zeros(1,length(N_values)));
Slc1=single(zeros(1,length(N_values)));
for i = 1:length(N_values)
    for j = 1:i
    Sl1(j)= single(x_val^(j-1)/factorial(j-1));
    end
    term1c = sort(Sl1,'ComparisonMethod', 'abs');
    Slc1(i)=single(sum(term1c));
end 

Sl11=single(zeros(1,length(N_values)));
Sld1=single(zeros(1,length(N_values)));

for i = 1:length(N_values)
    for j = 1:i
    Sl11(j)= single(x_val^(j-1)/factorial(j-1));
    end
    term2c= sort(sort(Sl11,'ComparisonMethod', 'abs'),'descend');
    Sld1(i)=single(sum(term2c));
end 



% Affichage des graphiques semilogy pour x=-5
figure;
semilogy(0:50, abs((Slc1 - somes(-5)))./abs(somes(-5)));
hold on;
semilogy(0:50, abs((Sld1 - somes(-5)))./abs(somes(-5)));


xlabel('N');
ylabel('Erreur absolue');
legend('Sc avec x=-5', 'Sd avec x=-5');
title('Erreur absolue de Sc(x;N) et Sd(x;N) par rapport à S(x;N)');


%Pour x=5

x_val2 = 5;
Sl2=single(zeros(1,length(N_values)));
Slc2=single(zeros(1,length(N_values)));
for i = 1:length(N_values)
    for j = 1:i
    Sl2(j)= single(x_val2^(j-1)/factorial(j-1));
    end
    term3c = sort(Sl2,'ComparisonMethod', 'abs');
    Slc2(i)=single(sum(term3c));
end 

Sl22=single(zeros(1,length(N_values)));
Sld2=single(zeros(1,length(N_values)));

for i = 1:length(N_values)
    for j = 1:i
    Sl22(j)= single(x_val2^(j-1)/factorial(j-1));
    end
    term4c= sort(sort(Sl22,'ComparisonMethod', 'abs'),'descend')
    Sld2(i)=single(sum(term4c));
end 

% Affichage des graphiques semilogy pour x=5
figure;
semilogy(0:50, abs((Slc2 - somes(5)))./abs(somes(5)));
hold on;
semilogy(0:50, abs((Sld2 - somes(5)))./abs(somes(5)));

xlabel('N');
ylabel('Erreur absolue');
legend('Sc avec x=5', 'Sd avec x=5');
title('Erreur absolue de Sc(x;N) et Sd(x;N) par rapport à S(x;N)');


function T=Erreur_relative(x)
T= zeros(1,51);
    for i=0:50
        T(i+1)=abs(exp(x)-sumes(x,i))/exp(x);
    end
end

function s=sumes(x,N)
    s=0;
    for i=0:N
        s=s+(x^i)/(factorial(i));
    end
end

function some= somes(x)
    some=zeros(1,51);
    for i=0:50
        some(i+1)=sumes(x,i);
    end
end

function borne=bornes(x)
borne=zeros(1,51);
    for N=0:50
            if x<0
                borne(N+1)= (exp(-x)*(-x).^(N+1))/factorial(N+1);
            else
                borne(N+1)= (x).^(N+1)/factorial(N+1);
            end
    end

end


