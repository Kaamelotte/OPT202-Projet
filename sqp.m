function [x,lme, lmi, info] = sqp(simul,x, lme, lmi, options)
################################################################################
## sqp :
## Renvoie les multiplicateurs de Lagrange et les abscisses et ordonnees des  
## noeuds
##
## INPUT  - simul : specification du simulateur
##        - x : vecteur contenant la valeur initiale
##        - lme : vecteur contenant les multiplcicateurs de Lagrange intiale
##        - options : structure speciafiant les parametres de fonctionnement 
##                    de l’algorithme
##
## OUTPUT - x : vecteur contenant la valeur finale
##        - lme : vecteur contenant les multiplicateurs de Lagrange pour les contrainte d egalite
##        - lmi : vecteur contenant les multiplicateurs de Lagrange pour les contrainte d inegalite
##
##        - info : structure donnant des informations sur le comportement du solveur
##
##############################################################################
	  
	##=== Verification de la consistance des arguments d'entree ========================##
	n = length(x);
	if mod(n,2) ~= 0
		info.status = 1 %length(x) impair
		return  
	end ##====================================================================##

	##=== Verification que l'utilisateur de demande pas plus d'info qu'il n'en veut=========##
	if options.rl == 1 && options.verb == 2
		options.verb = 1;
	end;##====================================================================##
	
	##=== Calcul des multiplicateurs lagrangiens d'egalite=============================##
    if length(lme) == 0
		[~,~,~,g,ae,~,~,indic] = simul(4,x,lme);
		lme = -ae'\g;
    end
	m = length(lme);
	##=======================================================================##

							#####################################
							##=== Algorithme de Newton ==========##
							#####################################
	nbIter = 0;
	nbSimul = 0;
	alpha = 1;  #pas initial = pas unite
	normFk = 1;
	#grdlp = 1;
  
	if (options.verb == 1)##=== Impression ===##
		fprintf('---------------------------------------------------------------------------------\n');
		fprintf('%4s %7s %10s %10s %10s %11s %9s %9s\n',...
				'iter','|gl|','|ce|','|x|','|lme|','alpha','phi','Q');
		fprintf('---------------------------------------------------------------------------------\n');
	end##================##  

##=== Boucle principale =======================================================##
	while true
		[~,ce,~,g,ae,~,~,indic] = simul(4,x,lme);

		nbSimul += 1;
		grdl = g + ae' * lme;
		
		##===Quotient des normes de Fk ============================================##
		normFkp = max(norm(ce,Inf),norm(grdl,Inf));
		Q = normFkp / normFk^2;
		normFk = normFkp;
		
		##=== Calcul de la direction de descente ===================================##
		[~,~,~,~,~,~,hl,~] = simul(5,x,lme);
		nbSimul += 1;
		dF = [ hl, ae'; ae, zeros(m,m) ];
		F = [ grdl ; ce ];
		
		dir = -dF\F; #(dk,muk)
		##===================================================================##
		
		#Faut lui trouver une petit place a ce petit chou
		phi = 0.5 * F' * F; #0.5*||F(z)||^2
		
		##===================================================================##
			
		if options.verb > 0 ##=== Impression ===##
			if options.verb == 2			
				fprintf('---------------------------------------------------------------------------------\n');
				fprintf('%4s %7s %10s %10s %10s %11s %9s %9s\n',...
						'iter','|gl|','|ce|','|x|','|lme|','alpha','phi','Q');
			end;
		  fprintf('%4d %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',...
				  nbIter, norm(grdl,inf),norm(ce,inf),norm(x,inf),norm(lme,inf),alpha,phi,Q);
		end ##================##    
			
							#####################################
							##=== Recherche lineaire ==============##
							#####################################
			
##=== Recherche lineaire pour le pas alpha ======================================##
		if options.rl == 0
			dphi = F' * dF;     #F(z)^T*F'(z)
			[alpha, nbSimul] = rl( x, lme, [], dir, simul, nbSimul, dphi, phi, options);
		end 
##=== Fin recherche lineaire ===================================================##
			
		##=== Calcul des nouveaux parametres pour Newton ==========================##
		x = x + alpha*dir(1:n);
		lme = lme + alpha*dir(n+1:length(dir));
			
		nbIter = nbIter + 1;
		info.niter = nbIter;
		##=====================================================================##

							#####################################
							##=== Test d arret====================##
							#####################################	
		
		##=== Test d'optimalite ===================================================##
		if (norm(grdl,inf) < options.tol(1)) && (norm(ce,inf) < options.tol(2))
			info.status = 0; #Solution trouvee
			info.niter = nbIter;
			break; #Sortie de Newton
		end##==================================================================##
			
		##===Test du nombre d'iterations deja effectuees==============================##
		if(nbIter >= options.maxit)
			info.status = 2; #Sortie de l'algo car pas de convergence
			info.tol(1) = norm(grdl,inf);
			info.tol(2) = norm(ce,inf);
            info.tol(3) = 0;
			break; 
		end ##==================================================================##
	end
##=== Fin Boucle principale ====================================================##
 
	if options.verb > 0##=== Impression ===##
		fprintf('---------------------------------------------------------------------------------\n');
	end ##================## 

	return
end  #Fin de la fonction
