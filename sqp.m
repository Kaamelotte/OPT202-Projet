function [x,lme,info] = sqp(simul,x, lme, lmi, options)
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
	m = length(lme);
	##=======================================================================##

							#####################################
							##=== Algorithme de Newton ==========##
							#####################################
	nbIter = 0;
	nbSimul = 0;
	alpha = 1;  #pas initial = pas unit�
	normFk = 1;
  
	if (options.verb == 1)##=== Impression ===##
		fprintf('---------------------------------------------------------------------------------\n');
		fprintf('%4s %7s %10s %10s %10s %11s %9s %9s\n',...
				'iter','|gl|','|ce|','|x|','|lme|','alpha','phi','Q');
		fprintf('---------------------------------------------------------------------------------\n');
	end##================##  

##=== Boucle principale =======================================================##
	while true		
		[e,ce,ci,g,ae,ai,~,indic] = simul(4,x,lme);
		nbSimul += 1;
		grdl = g + ae' * lme;
		
		##===Quotient des normes des derives de Fk==================================##
		normFkp = max(norm(ce,Inf),norm(grdl,Inf));
		Q = normFkp / normFk^2;
		normFk = normFkp;
		
		##=== Calcul de la direction de descente ===================================##
		[~,~,~,~,~,~,hl,~] = simul(5,x,lme, []);
		nbSimul += 1;
		dF = [ hl, ae'; ae, zeros(m,m) ];
		F = [ g ; ce ];
		dir = -dF\F;    
		##===================================================================##
			
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
							##=== Recherche lin�aire ==============##
							#####################################
							
##=== Recherche lin�aire pour le pas alpha ======================================##
		if options.rl == 0
			dphi = F' * dF;     #F(z)^T*F'(z)
			[alpha, nbSimul] = rl(simul,x, nbSimul, lme, dir, dphi, phi, 1e-4, options);
		end 
##=== Fin recherche lin�aire =============================================##
			
		##=== Calcul des nouveaux param�tres pour Newton ========================##
		x = x + alpha*dir(1:n);
		lme = (1-alpha)*lme + alpha*dir(n+1:length(dir));
			
		nbIter = nbIter + 1;
		info.niter = nbIter;
		##===================================================================##
    
    if contraintes d'in�galit�
      pb osculateur
    end
		
							#####################################
							##=== Test d arret====================##
							#####################################	
		
		##=== Test d'optimalit� =================================================##
		if (norm(grdl,inf) < options.tol(1)) && (norm(ce,inf) < options.tol(2))
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
 
	if options.verb > 0##=== Impression ===##
		fprintf('---------------------------------------------------------------------------------\n');
	end ##================## 

	return
end  #Fin de la fonction
