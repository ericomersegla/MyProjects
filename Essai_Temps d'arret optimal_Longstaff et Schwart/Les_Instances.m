function N=Les_Instances()
    ins=input('Entrer le nombre pour la periode ');
    if ins<=18 %soit 18 mois ,nous considerons les instances comme le nombre de mois
        display('Petite instance');
        N=ins;
    elseif 18<ins& ins<=72 %soit 1 an 6 mois à 6 ans 
        display('Moyenne instance');
        N=ins;
    else  
        display('Grande instance'); %soit plus de 6 ans 
        N=ins;
    end
      