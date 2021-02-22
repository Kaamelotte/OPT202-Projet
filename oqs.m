function [x,lme, lmi, info] = oqs(simul,x, lme, lmi, options)
################################################################################
## sqp :                                                                      ##
## Renvoie les multiplicateurs de Lagrange et les abscisses et ordonn�es des  ##
## noeuds                                                                     ##
##                                                                            ##
## INPUT  - simul : sp�cification du simulateur                               ##
##        - x : vecteur contenant la valeur initiale                          ##
##        - lme : vecteur contenant les multiplcicateurs de Lagrange intiale   ##
##        - options : structure sp�ciafiant les param�tres de fonctionnement  ##
##                    de l�algorithme                                         ##
##                                                                            ##
## OUTPUT - x : vecteur contenant la valeur finale                            ##
##        - lme : vecteur contenant les multiplcicateurs de Lagrange finaux    ##                                                      ##
##        - info : structure donnant des informations sur le comportement du  ##
##                  solveur                                                   ##
##                                                                            ##
##############################################################################
	  
	##=== V�rification de la consistance des arguments d'entr�e ========================##
	n = length(x);
	if mod(n,2) ~= 0
		info.status = 1 %length(x) impair
		return  
	end ##====================================================================##

	##=== Calcul des multiplicateurs lagrangiens d'�galit�=============================##
    if length(lme) == 0 || length(lmi) == 0
		[~,~,~,g,ae,ai,~,indic] = chs(4,x,lme, lmi);
		lme = -ae'\g;
		lmi = -ai'\g;
    end
	me = length(lme);
	mi = length(lmi);
	##=======================================================================##

							#####################################
							##=== Algorithme de Newton ==========##
							#####################################
	nbIter = 0;
	nbSimul = 0;

##=== Boucle principale =======================================================##
	while true
		[~,ce,ci,g,ae,ai,~,indic] = simul(4,x,lme);
		nbSimul += 1;
		
		##=== Calcul de la direction de descente ===================================##
		[~,~,~,~,~,~,hl,~] = simul(5,x,lme, []);
		nbSimul += 1;
		 [L, d, flag] = cholmod(hl, 1.e-5, 1.e+5);
		 M = L*diag(d)*L';
		 		 		 
		 #[d, obj, information, lm] = qp (X0, H, Q,  A,    B, LB, UB, A_LB, A_IN, A_UB)
		 [d, obj, information, lm] = qp ( d, M,  g, ae, -ce, [] ,  [] ,     []  ,    ai  ,   -ci   );
		
		##=====================================================================##
			
		##=== Calcul des nouveaux param�tres pour Newton ========================##
		x = x + d;
		lme = lm(1:me); #ordre ?
		lmi = lm(me+1:me +mi);
		
		nbIter = nbIter + 1;
		info.niter = nbIter;
		##===================================================================##

							#####################################
							##=== Test d arret====================##
							#####################################	
		
		##=== Test d'optimalit� =================================================##		
		[~,ce,ci,g,ae,~,~,indic] = simul(4,x,lme);
		nbSimul += 1;
		grdl = g + ae' * lme;
		
		if (norm(grdl,inf) < options.tol(1) && norm(ce,inf) < options.tol(2) && norm(lmi-ci,inf) < option.tol(3) )
			info.status = 0; #Solution trouv�e
			info.niter = nbIter;
			break; #Sortie de Newton
		end##================================================================##
			
		##===Test du nombre d'it�rations d�j� effectu�es============================##
		if(nbIter > options.maxit)
			info.status = 2; #Sortie de l'algo car pas de convergence
			break; 
		end ##================================================================##
	end
##=== Fin Boucle principale ====================================================##
 
	return
end  #Fin de la fonction
