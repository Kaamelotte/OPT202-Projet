function [] = checker(castest, test = 1)
	global A B L R S color lgd;
	if test != "-1"
		[L,xy,A,B,R,S] = casTest(castest)
	end
	switch test
		case '1' ##=== Trace de la chaine ===##	  
			disp('trace de la chaine');
			hold on;
			chs(1,xy);
			chs(6,xy);
			hold off;
		case 'c'  ##=== Energie potentielle et contraintes ======##
			disp("calcul de l'energie et des contraintes");
			[e,ce,ci,~,~,~,~,indic] = chs(2,xy)
		case 'g'  ##=== Gradient de e et jacobienne de c  =================##
			disp("calcul du gradient de e et de la jacobienne de c");
			[~,~,~,g,ae,ai,~,indic] = chs(4,xy);
			g
			ae = full(ae)
			ai = full(ai)
		case 'hl' ##=== Hessien du lagrangien ============##
			disp("calcul du hessien du lagrangien");
			[~,~,~,~,~,~,hl,indic] = chs(5,xy);
			hl = full(hl)
			valuerPropre = eig(hl)
			sum(eig(hl)>= 0)/length(eig(hl)) 
		case 'grad' ##=== Verification du gradient de l'energie potentielle ===##
			disp("verification du gradient de l'energie potentielle"); 
			verifierGradient(xy);
		case 'lm' ##=== Multiplicateur de lagrange===##
			[~,~,~,g,ae,ai,~,indic] = chs(4,x);
			lme = -ae'\g
			lmi = -ai'\g
		case 'L' ##=== longueur de la chaine ===##
			[~,ce,~,~,~,~,~,indic] = chs(2,xy);
			longueur = L
			longueur_actuelle = sqrt( ce + L.^2 )
		case 'cholmeod' ##=== Test de la fonction cholmeod ===##
			disp("Test de la fonction cholmeod\n"); 
			#Tolerance
			tol_small = 1.e-5; 
			tol_big = 1.e+5;
			[~,~,~,~,~,~,hl,indic] = chs(5,xy);
			[L, d, flag] = cholmeod(hl, tol_small, tol_big);
			full(hl)
			full(L)
			D = diag(d)
			M = L*D*L'
			disp("M est bien defini positive\n")
			eig(M) #valeur propre de M
			disp("Erreur d'approximation de hl\n")
			E =  norm(hl - M,inf)
	end##============================================================================##
	
	return
end