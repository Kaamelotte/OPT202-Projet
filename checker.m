function [] = checker(castest, test = 1)
	global A B L R S color lgd;
	if strcmp(test, "-1")
		return
	end
	[L,xy,A,B,R,S] = casTest(castest)
	switch test
		case '1'   
			disp('trace de la chaine');
			hold on;
            chs(6,xy);
			#chs(1,xy);
			hold off;
		case 'c'
			disp("calcul de l'energie et des contraintes");
			[e,ce,ci,~,~,~,~,indic] = chs(2,xy)
		case 'g' 
			disp("calcul du gradient de e et de la jacobienne de c");
			[~,~,~,g,ae,ai,~,indic] = chs(4,xy);
			g
			ae = full(ae)
			ai = full(ai)
    case 'grdl'
        size_xy = size(xy)
            [~,~,~,g,ae,ai,~,indic] = chs(4,xy);
            size_g = size(g)
            size_ae = size(ae)
            size_ai = size(ai)
		case 'hl' 	
			disp("calcul du hessien du lagrangien");
			[~,~,~,~,~,~,hl,indic] = chs(5,xy);
			hl = full(hl)
			valeurPropre = eig(hl)
			sum(eig(hl)>= 0)/length(eig(hl)) 
		case 'grad'
			disp("verification du gradient de l'energie potentielle"); 
			verifierGradient(xy);
        case 'lm' 
            disp("Multiplicateurs de lagrange")
			[~,~,~,g,ae,ai,~,indic] = chs(4,x);
			lme = -ae'\g
			lmi = -ai'\g
        case 'L'
            disp("Longueur de la chaine")
			[~,ce,~,~,~,~,~,indic] = chs(2,xy);
			longueur = L
			longueur_actuelle = sqrt( ce + L.^2 )
		case 'cholmod' 
			disp("Test de la fonction cholmod\n"); 
			#Tolerance
			tol_small = 1.e-5; 
			tol_big = 1.e+5;
			[~,~,~,~,~,~,hl,indic] = chs(5,xy,[],[]);
			[L, d, flag] = cholmod(hl, tol_small, tol_big);
			disp('hl=')
            disp(full(hl))
            disp('L=')
			disp(full(L))
			D = diag(d)
			M = L*D*L'
			disp("M est bien definie positive\nSes valeurs propres sont:")
			disp(eig(M)) #valeur propre de M
			disp("Erreur d'approximation de hl:\n")
			E =  norm(hl - M,inf)
	end##====================================================================##
	
	return
end