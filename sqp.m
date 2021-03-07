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

	##=== Calcul des multiplicateurs lagrangiens ====================================##
	[~,ce,ci,g,ae,ai,~,indic] = chs(4,x,lme, lmi);
    if length(ai) > 0 && length(lmi) == 0
        me = size(ae)(1);       mi = size(ai)(1);         m = size([ae; ai])(1);
        lm0 = zeros(m,1);       H = [ae;ai]*[ae;ai]';   q = [ae;ai]*g;
        A = [zeros(1,me), ci']; b = 0;
        lb = zeros(m,1);          ub = repmat(Inf,m,1);
        A_lb = ci;                    A_in = zeros(mi,m);   A_ub = repmat(Inf,mi,1);
	    [lm, obj, information, l] = qp (lm0, H, q, A, b, lb, ub, A_lb, A_in, A_ub);
        #lm = -[ae; ai]' \ g
        lme =  lm(1:me);
        lmi = lm(me+1:m);
    elseif length(lme) == 0
	    lme = -ae'\g;
    end
	me = length(lme);
	mi = length(lmi);
	##=======================================================================##

	info.niter = 0;
	info.nbSimul = 0;
	alpha = 1;  #pas initial = pas unite
	xm = x;
  
	if (options.verb == 1) ##=== Impression ===##
		fprintf('%4s %7s %10s %10s %10s %10s %10s %11s %10s',...
				'iter','|gl|','|ce|','|ci|','|x|','|lme|','|lmi|'); fprintf('\n');
	end ##==============================##  

##=== Boucle principale =======================================================##
	while true
		[e,ce,ci,g,ae,ai,~,indic] = simul(4,x,lme, lmi);		
		[~,~,~,~,~,~,hl,~] = simul(5,x,lme, lmi);
		info.nbSimul += 2;
		grdl = g + [ae ; ai ]' * [lme ; lmi ];
		#grdl = g + ae'*lme 
		
        ##=== Test d'optimalite ====================================================##
		if  (norm(ce,inf) < options.tol(2))  && (norm(min(lmi,-ci),inf) < options.tol(3) ) #&&(norm(grdl,inf) < options.tol(1)) 
			info.status = 0; #Solution trouvee
			break; #Sortie de la boucle principale
		end##==================================================================##
		
		##=== Calcul de la fonction F et de sa derive ===================================##
		if options.quad == 0 || options.rl == 0 ;
			dF = [ hl, [ae;ai]' ; [ae; ai], zeros(me+mi,me+mi) ];
			F = [ grdl ; [ce; ci] ];
			FPQ = [ g ; [ce ; ci ] ];
		end ##=================================================================##
		
        ##=== Calcul de la direction de descente ======================================##
		if options.quad == 0          ##=== Algorithme de Newton ==========##
			dirct = -dF\FPQ; #(dk,lm_PQ)
		elseif options.quad == 1    ##=== Algorithme de Josephy-Newton ===##
			if options.deriv == 2
				[L, d, flag] = cholmod(hl, 1.e-5, 1.e+5);
				M = L*diag(d)*L';
			elseif (info.niter == 0) &&(options.deriv == 1)
				M = eye(n); # matrice initiale
			end
			#[d, obj, information, lm] = qp (X0            , H, Q,  A,    B, LB, UB, A_LB, A_IN, A_UB)
			[dk, obj, information, lm] = qp (zeros(n,1), M,  g, ae, -ce, [] ,  [] ,     []  ,    ai  ,   -ci   );
			if information.info != 0
				info.status = information.info
				break; #on sort de la boucle while principale
			end
			dirct = [dk; lm ];
		end
        dk = dirct(1:n);
		lmePQ = dirct(n+1:n+me);
		lmiPQ = dirct(n+me+1:length(dirct));
        ##=== Fin calcul de la dirctection de descente =================================##
        
        if (options.verb = 4) && (mod(info.niter,5) == 0) ##=== Plot ===##
            figure(2)
            chs(6,x,lme,lmi);
            chs(1,x,lme,lmi);
            hold on;
            quiver(x(1:n/2),x(n/2+1:n),dk(1:n/2),dk(n/2+1:n));
            title(["itération ", num2str(info.niter)]);
            pause();
            hold off;
        end
        
		if options.verb > 0 ##=== Impression ===##
			if options.verb == 2			
				fprintf('-------------------------------------------------------------------------------------------\n');
				fprintf('%4s %7s %10s %10s %10s %10s %10s\n',...
						'iter','|gl|','|ce|','|ci|','|x|','|lme|','|lmi|');
			end;
		    fprintf('%4d %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',...
			    info.niter, norm(grdl,inf),norm(ce,inf),norm(ci,inf),norm(x,inf),norm(lme,inf),norm(lmi,inf));
		end
			
        
		if options.rl == 0 ##=== Recherche lineaire pour le pas alpha ===##
			phi = 0.5 * F' * F; #0.5*||F(z)||^2
			dphi = F' * dF;     #F(z)^T*F'(z)
			[alpha, info] = rl( x, lme, lmi, dirct, simul, info, dphi, phi, options);
        elseif info.niter > 50
            alpha = 0.5;
		end 
		
		##=== Nouveaux parametres pour Newton ===================================##
		xm = x;
        x = x + alpha*dk;
		lme = lme +  alpha*(lmePQ - lme);
		lmi = lmi +  alpha*(lmiPQ -lmi);
		
		if (options.quad == 1) && (options.deriv == 1) ### === BFGS =====================##
			[M, info] = bfgs(simul, info, M, xm, lme, lmi, dirct, grdl, options );
		end 
        
		##=== Test du nombre d'iterations deja effectuees =============================##
		info.niter = info.niter + 1;
        if(info.niter >= options.maxit)
			info.status = 2; #Sortie de l'algo car pas de convergence
			info.tol(1) = norm(grdl,inf);
			info.tol(2) = norm(ce,inf);
            info.tol(3) = norm(min(lmi,-ci),inf);
			break; 
		end ##=================================================================##
	end ##=== Fin Boucle principale ==============================================##
    #close(2);
	return
end  #Fin de la fonction
