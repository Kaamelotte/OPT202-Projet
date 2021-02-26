function [x,lme, lmi, info] = sqp(simul,x, lme, lmi, options)
#############################################################################
## sqp :
## Renvoie les multiplicateurs de Lagrange et les abscisses et ordonnees des  
## noeuds
##
## INPUT  
##        - simul : specification du simulateur
##        - x : vecteur contenant la valeur initiale
##        - lme : vecteur contenant les mult de Lagrange pour les contraintes d'égalité
##        - lmi : vecteur contenant les mult de Lagrange pour les contraintes d'inégalité
##        - options : structure specifiant les parametres de fonctionnement de l’algorithme
## OUTPUT 
##        - x : vecteur contenant la valeur finale
##        - lme : vecteur contenant les multiplicateurs de Lagrange pour les contrainte d egalite
##        - lmi : vecteur contenant les multiplicateurs de Lagrange pour les contrainte d inegalite
##        - info : structure donnant des informations sur le comportement du solveur
##
#############################################################################
	  
	##=== Verification de la consistance des arguments d'entree =======================##
	n = length(x);
	if mod(n,2) != 0
		info.status = 1 %length(x) impair
		return  
	end ##===================================================================##

	##=== Verification que l'utilisateur de demande pas plus d'info qu'il n'en veut ==========##
	if options.rl == 1 && options.verb == 2
		options.verb = 1;
	end; ##===================================================================##
	
	##=== Calcul des multiplicateurs lagrangiens ====================================##
	if length(lme) == 0 
		    [~,~,~,g,ae,ai,~,indic] = chs(4,x,lme, lmi);
		    lme = -ae'\g;
	end
	if length(lmi) == 0
		[~,~,~,g,ae,ai,~,indic] = chs(4,x,lme, lmi);
		lmi = -ai'\g;
	end
	me = length(lme);
	mi = length(lmi);
	##=======================================================================##

	info.niter = 0;
	nbSimul = 0;
	alpha = 1;  #pas initial = pas unite
	normFk = 1;
	#grdlp = 1;
  
	if (options.verb == 1) ##=== Impression ===##
		fprintf('---------------------------------------------------------------------------------\n');
		fprintf('%4s %7s %10s %10s %10s %10s %10s %11s %9s\n',...
				'iter','|gl|','|ce|','|ci|','|x|','|lme|','|lmi|','alpha','Q');
		fprintf('---------------------------------------------------------------------------------\n');
	end ##==============================##  

##=== Boucle principale =======================================================##
	while true
		#=== Mise a jour de donnees ===============================================##
		[~,ce,ci,g,ae,ai,~,indic] = simul(4,x,lme, lmi);		
		[~,~,~,~,~,~,hl,~] = simul(5,x,lme, lmi);
		nbSimul += 2;
		
		grdl = g + [ae ; ai ]' * [lme ; lmi ];
		
		##=== Test d'optimalite ===================================================##
		if (norm(grdl,inf) < options.tol(1)) && (norm(ce,inf) < options.tol(2))  && (norm(min(lmi,-ci),inf) < options.tol(3) )
			info.status = 0; #Solution trouvee
			break; #Sortie de la boucle principale
		end##==================================================================##
		
		##=== Quotient des normes de Fk ===========================================##
		normFkp = max(norm(ce,Inf),norm(grdl,Inf));
		Q = normFkp / normFk^2;
		normFk = normFkp;
		
		##=== Calcul de la fonction F et de sa derive==================================##
		if options.quad == 0 || options.rl == 0 ;
			dF = [ hl, [ae;ai]' ; [ae; ai], zeros(me+mi,me+mi) ];
			F = [ grdl ; [ce; ci] ];
			Fpq = [ g ; [ce ; ci ] ];
			F = Fpq;
		end ##=================================================================##
		
##=== Calcul de la direction de descente =========================================##
		if options.quad == 0 			##=== Algorithme de Newton ==========##
			dir = -dF\F; #(dk,muk)
			
		elseif options.quad == 1		##=== Algorithme de Josephy-Newton nom ????? ===##
			 [L, d, flag] = cholmod(hl, 1.e-5, 1.e+5);
			 M = L*diag(d)*L';
			 #[d, obj, information, lm] = qp (X0             , H, Q,  A,    B, LB, UB, A_LB, A_IN, A_UB)
			 [dk, obj, information, lm] = qp (ones(n,1) , M,  g, ae, -ce, [] ,  [] ,     []  ,    ai  ,   -ci   );
			if information.info != 0
				info.status = information.info;
				break;
			end
			dir = [dk; lm ];
		end
##=== Fin calcul de la direction de descente =======================================##		
			
		if options.verb > 0 ##=== Impression ===##
			if options.verb == 2			
				fprintf('---------------------------------------------------------------------------------\n');
				fprintf('%4s %7s %10s %10s %10s %11s %9s\n',...
						'iter','|gl|','|ce|','|x|','|lme|','alpha','Q');
			end;
		  fprintf('%4d %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e  \n',...
				  info.niter, norm(grdl,inf),norm(ce,inf),norm(x,inf),norm(lme,inf),alpha,Q);
		end ##============================##    
			
##=== Recherche lineaire pour le pas alpha ========================================##
		if options.rl == 0
			phi = 0.5 * F' * F; #0.5*||F(z)||^2
			dphi = F' * dF;     #F(z)^T*F'(z)
			[alpha, nbSimul] = rl( x, lme, lmi, dir, simul, nbSimul, dphi, phi, options);
		end 
##=== Fin recherche lineaire ====================================================##
		
		##=== Calcul des nouveaux parametres pour Newton ============================##
		x = x + alpha*dir(1:n);
		lme = lme + alpha*(dir(n+1:n+me)-lme);
		lmi = lmi + alpha*(dir(n+me+1:length(dir))-lmi);
		
		info.niter = info.niter + 1;
		##=====================================================================##

		##=== Test du nombre d'iterations deja effectuees ==============================##
		if(info.niter >= options.maxit)
			info.status = 2; #Sortie de l'algo car pas de convergence
			info.tol(1) = norm(grdl,inf);
			info.tol(2) = norm(ce,inf);
            info.tol(3) = norm(min(lmi,-ci),inf);
			break; 
		end ##=================================================================##
	end
##=== Fin Boucle principale ====================================================##
 
	if options.verb > 0 ##=== Impression ===##
		fprintf('---------------------------------------------------------------------------------\n');
	end ##============================## 

	return
end  #Fin de la fonction
