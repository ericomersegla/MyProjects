/* Importation de la base de données et création de la table Simul_survie */

PROC IMPORT DATAFILE="/home/u64088625/sasuser.v94/CREDIT_RISK/simul_survie.xlsx"
    OUT=Simul_survie
    DBMS=XLSX
    REPLACE;
    GETNAMES=YES;
RUN;


/* Utilisation de la table imul_survie */
DATA Analyse_survie;
    SET Simul_survie;
RUN;

/*Affichage des 10 premières lignes*/
PROC PRINT DATA=Analyse_survie (OBS=10);
RUN;

/*ID_compte	date_observation	date_defaut	date_fermeture	age	revenu	Jrs_delinquance	chomage*/
/*Les stats descritives*/
PROC UNIVARIATE DATA=Analyse_survie;
    VAR age	revenu	Jrs_delinquance	chomage;/* Analyse toutes les variables numériques */
    HISTOGRAM / NORMAL;  /* Affiche aussi un histogramme avec la loi normale */
RUN;

/*Les stats descritives 2*/
PROC MEANS DATA=Analyse_survie N MIN MEAN MEDIAN MAX ;
    VAR age	revenu	Jrs_delinquance	chomage;  /* Analyse toutes les variables numériques */
RUN;

/* Traitement de données pour la modelisation du risque credit*/

%let date_fin="31DEC2023"d;
data Analyse_survie;
	 SET Simul_survie;
		/*"Censure": qui indique l'occurence du defaut sur la durée de vie*/
		Censure=(not missing(date_defaut));
		
		/*le temps de survie*/
		if not missing(date_defaut) then
			T_survie=intck("months",date_observation,date_defaut);
		else if not missing(date_fermeture) then
			T_survie=intck("months",date_observation,min(date_fermeture,&date_fin.));
		else T_survie=intck("months",date_observation,&date_fin.);
		
		/*Suppression de temps de survie nul*/
		if T_survie=0 then 
				delete;
				
		/* Linearisation du revenu*/
		L_revenu=log(revenu);
run;
		
	/* REGRESSION AVEC PROC LIFEREG*/
proc lifereg data=analyse_survie;
	model T_survie*censure(0) = age	Jrs_delinquance	chomage L_revenu;
run;
		
		/* REGRESSION AVEC PROC PHEREG POUR LA SÉLECTION DES VARIABLES*/
proc phreg data=analyse_survie;
	model T_survie*censure(0) = age	Jrs_delinquance	chomage L_revenu/ selection=backward;
run;	
		

	/* REGRESSION AVEC PROC LIFEREG SANS LA VARIABLE AGE--- AVEC LA LOI EXPONENTIELLE*/
proc lifereg data=analyse_survie;
	model T_survie*censure(0) =	chomage Jrs_delinquance L_revenu/DISTRIBUTION=exponential;
run;


	/* REGRESSION AVEC PROC LIFEREG SANS LA VARIABLE AGE--- AVEC LA LOI loglogistique*/
proc lifereg data=analyse_survie;
	model T_survie*censure(0) =	chomage Jrs_delinquance L_revenu/distribution=llogistic;
run;



	/* REGRESSION AVEC PROC LIFEREG SANS LA VARIABLE AGE--- AVEC LA LOI lognormale*/
proc lifereg data=analyse_survie;
	model T_survie*censure(0) =	chomage Jrs_delinquance L_revenu/DISTRIBUTION=lnormal;
run;


	/* REGRESSION AVEC PROC LIFEREG SANS LA VARIABLE AGE--- AVEC LA LOI weibull*/
proc lifereg data=analyse_survie;
	model T_survie*censure(0) =	chomage Jrs_delinquance L_revenu/distribution=weibull;
run;


/* MODELISATION AVEC LOGISTIC*/

proc logistic data=analyse_survie descending;
   model censure(event='0') = age chomage Jrs_delinquance L_revenu;
run;










