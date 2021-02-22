function [x,lme, lmi, info] = oqs(simul,x, lme, lmi, options)

	##=== Vérification de la consistance des arguments d'entrée ========================##
	n = length(x);
	if mod(n,2) ~= 0
		info.status = 1 %length(x) impair
		return  
	end ##====================================================================##
	
	##=== Verification que l'utilisateur de demande pas plus d'info qu'il n'en veut=========##
	if options.verb == 2
		options.verb = 1;
	end;##====================================================================##
	

	##=== Calcul des multiplicateurs lagrangiens d'égalité=============================##
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
	
	if (options.verb == 1)##=== Impression ===##
		fprintf('---------------------------------------------------------------------------------\n');
		fprintf('%4d %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e \n',...
					'iter','|ce|','|ci|','|x|','|lme|','|lmi-ci|', '|E|');
		fprintf('---------------------------------------------------------------------------------\n');
	end##================##  
	
##=== Boucle principale =======================================================##
	while true
		[~,ce,ci,g,ae,ai,~,indic] = simul(4,x,lme);
		nbSimul += 1;
		
		##=== Calcul de la direction de descente ===================================##
		[~,~,~,~,~,~,hl,~] = simul(5,x,lme, []);
		nbSimul += 1;
		 [L, d, flag] = cholmod(hl, 1.e-5, 1.e+5);
		 M = L*diag(d)*L';
		grdl = g + ae' * lme;
		 		 		 
		 #[d, obj, information, lm] = qp (X0, H, Q,  A,    B, LB, UB, A_LB, A_IN, A_UB)
		 [d, obj, information, lm] = qp ( d, M,  g, ae, -ce, [] ,  [] ,     []  ,    ai  ,   -ci   );
		
		##=====================================================================##
			
		##=== Calcul des nouveaux paramètres pour Newton ========================##
		x = x + d;
		lme = lm(1:me); #ordre ?
		lmi = lm(me+1:me +mi);
		
		nbIter = nbIter + 1;
		info.niter = nbIter;
		##===================================================================##

		[~,ce,ci,g,ae,~,~,indic] = simul(4,x,lme);
		nbSimul += 1;
		
		if options.verb > 0 ##=== Impression ===##
		  fprintf('%4d %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e \n',...
				  nbIter, norm(ce,inf),norm(ci,inf),norm(x,inf),norm(lme,inf),norm(lmi,inf), norm(hl - M,inf));
		end ##================##    
		
							#####################################
							##=== Test d arret====================##
							#####################################	
		
		##=== Test d'optimalité =================================================##		
		if (norm(grdl,inf) < options.tol(1) && norm(ce,inf) < options.tol(2) && norm(lmi-ci,inf) < option.tol(3) )
			info.status = 0; #Solution trouvée
			info.niter = nbIter;
			break; #Sortie de Newton
		end##================================================================##
			
		##===Test du nombre d'itérations déjà effectuées============================##
		if(nbIter > options.maxit)
			info.status = 2; #Sortie de l'algo car pas de convergence
			break; 
		end ##================================================================##
	end
##=== Fin Boucle principale ====================================================##
 
	return
end  #Fin de la fonction
