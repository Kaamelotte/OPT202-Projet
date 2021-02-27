function [M, info] = bfgs(simul, info, M, xm, lme, lmi, dir, alpha, grdl )
	n = length(xm);
	#variation de x
	delta = alpha * dir(1:n);
	#variation du gradient
	[~,~,~,gm,aem,aim,~,indic] = simul(4,xm,lme, lmi);
	info.nbSimul += 1;				
	grdlbis = gm + [aem ; aim ]' * [lme ; lmi ];

	gammaL = grdl - grdlbis;

	#theta
	theta = 1;
	D = delta' * M * delta;
	fprintf( " delta'*M*delta : %e\n", D);
	S = gammaL' * delta;
	if S < 0.2 * D
		theta = 0.8 * D / (D - S);
	end
	#correction de powel
	gammaK = (1-theta) * M * delta + theta * gammaL;

	eta = (gammaK'*gammaK) / (gammaK' * delta);
	if info.niter == 0 
		M = eta*eye(n);
	end
	M = M - (M*delta*delta' * M )/ D  + eta;
	return
end