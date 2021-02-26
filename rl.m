function [alpha, nbSimul] = rl(x, lme, lmi, dir, simul, nbSimul, dphi, phi, options, nMax = 10, omega = 1e-4)
#############################################################################
## rl :
## Renvoie le alpha optimal trouve par recherche lineaire
##
## INPUT  
##        - x : vecteur contenant la valeur courante a optimiser
##        - lme : vecteur contenant les mult de Lagrange courants pour les contraintes d'egalite
##        - lmi : vecteur contenant les mult de Lagrange courants pour les contraintes d'inegalite
##        - dir = (dk,muk) = dir de descente du gradient pour x et lambda
##        - simul : specification du simulateur
##        - nbSimul: nombre d'appel du simulateur avant la recherche lineaire
##        - dphi: derivee de phi 
##        - phi: fonction de recherche lineaire
##        - options : structure specifiant les parametres de fonctionnement de l algorithme
## OUTPUT 
##        - alpha: alpha optimal trouve
##        - nbSimul: nbre d'appel au simulateur apres la recherche lineaire
#############################################################################

	n = length(x);
	me = length(lme);
	mi = length(lmi);
	
    dk = dir(1:n);#directionn de descente pour x
	lmPQ = dir(n+1:length(dir)); 
	muk = lmPQ - [lme ; lmi ]; # direction de descente pour les lambda
	
	i = 0; #nb d'iterations de la RL
	alpha = 1;
	pente = omega * dot(dphi,[dk;muk]); #omega.phi'(z_k).p_k	
	
	
	
	if options.verb == 2 ##=== Impression ===##
		fprintf('---------------------------------------------------------------------------------\n');
		fprintf("  Recherche lineaire d'Armijo: |d| = %.2e\n",norm(dir,inf));
		fprintf('    Simul =%d, phi = %.5e, omega*pente = %.5e\n\n',nbSimul, phi,omega*pente);
		fprintf('    %10s %15s %13s\n','alpha','phip-phi','DF(phi)');
	end ##============================## 
	
##=== Boucle sur i ============================================================##
	while true
		xp = x + alpha*dk; 			
		lmep = lme + alpha*muk(1:me); #lambda_{k+1} = \lambda_k + \alpha * mu_k)
		lmip =  lmi  + alpha*muk(me+1:me +mi);
		
		[~,cep,cip,gp,aep,aip,~,~] = simul(4,xp,lmep, lmip);
		nbSimul += 1;
		
		grdlp = gp + [aep ; aip ]' * [lmep ; lmip ];
		
		Fp = [ grdlp ; [cep; cip] ]; # [ grad l ; c ] cf 2.3 TP2
        phip = 0.5 * Fp' * Fp; #0.5*||F(z_k + alpha_k.p_k)||^2}
		
		##=== Test d'optimalite sur phi =======================================##
		if phip <= phi + alpha*pente # (phip-phi)/alpha <= pente
			break; #Sortie de Recherche lineaire
		end ##===========================================================##
		
		if options.verb == 2 ##=== Impression ===##
			fprintf('    %12.4e %14.5e %14.5e\n',alpha,phip-phi,(phip-phi)/alpha);
		end ##============================## 
		
		alpha = alpha/2;
		i = i+1;
		##=== Test du nombre d'iterations deja effectuees ==========================##
		if(i > nMax)
			info.status = 3; #PB dans la RL
			break; 
		end ##==============================================================##
	end
##=== Fin boucle i ============================================================##

	if (i!= 0) && (options.verb == 2) ##=== Impression ===##
		fprintf('  |gl| = %.3e, |ce| = %.3e\n',norm(grdlp,inf),norm(cep,inf));
	end ##============================## 
	
	return
end