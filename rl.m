function [alpha, nbSimul] = rl(x, lme, lmi, dir, simul, nbSimul, dphi, phi, options, omega = 1e-4)#dir = (dk,muk) = dir de descente du gradient pour x et lambda	n = length(x);	me = length(lme);	mi = length(lmi);		i = 0; #nb d'it�rations de la RL	alpha = 1;	pente = omega * dot(dphi,dir); #omega.alpha_k.phi'(z_k).p_k			lm = dir(n+1:length(dir)); # descente gradient pour les lambda	dir = dir(1:n);		if options.verb == 2 ##=== Impression ===##		fprintf('---------------------------------------------------------------------------------\n');		fprintf("  Recherche lin�aire d'Armijo: |d| = %.2e\n",norm(dir,inf));		fprintf('    Simul =%d, phi = %.5e, pente = %.5e\n\n',nbSimul, phi,pente);		fprintf('    %10s %15s %13s\n','alpha','phip-phi','DF(phi)');	end ##================## 	##=== Boucle sur i ============================================================##	while true		xp = x + alpha*dir; 		#lp = 		#lmp = (1-alpha)*lm + alpha*dir(n+1:length(dir));							lmep = lme + alpha* lm(1:me); #lambda_{k+1} = \lambda_k + \alpha * \mu_k		lmip =  lmi  + alpha* lm(me+1:me +mi);				[~,cep,cip,gp,aep,aip,~,~] = simul(4,xp,lmep, lmip);		nbSimul += 1;				grdlp = gp + [aep ; aip ]' * [lmep ; lmip ];				Fp = [ grdlp ; [cep; cip] ]; # [ grad f ; c ] cf 2.3 TP2		phip = 0.5 * Fp' * Fp; #0.5*||F(z_k+alpha_k.p_k)||^2		if options.verb == 2##=== Impression ===##			fprintf('    %12.4e %14.5e %14.5e\n',alpha,phip-phi,(phip-phi)/alpha);		end##================## 				##=== Test d'optimalit� sur phi =======================================##		if phip <= phi + alpha*pente # (phip-phi)/alpha <= pente			break; #Sortie de Recherche lin�aire		end##============================================================##				alpha = alpha/2;		i = i+1;		##===Test du nombre d'it�rations d�j� effectu�es==========================##		if(i > 10)			info.status = 3; #PB dans la RL			break; 		end##==============================================================##	end##=== Fin boucle i ============================================================##	if options.verb == 2##=== Impression ===##		fprintf('  |gl| = %.3e, |ce| = %.3e\n',norm(grdlp,inf),norm(cep,inf));	end##================## 		returnend;