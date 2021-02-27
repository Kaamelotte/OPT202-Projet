function [M, info] = bfgs(simul, info, M, xm, lme, lmi, dir, grdl, options )
	
	n = length(xm);
	#variation de x
	delta = dir(1:n);
	#variation du gradient
	[~,~,~,gm,aem,aim,~,indic] = simul(4,xm,lme, lmi);
	info.nbSimul += 1;				
	grdlbis = gm + [aem ; aim ]' * [lme ; lmi ];

	gammaL = grdl - grdlbis;

	#theta
	theta = 1;
	D = delta' * M * delta;
	S = gammaL' * delta;
	Q = D / (D - S);
	if S < 0.2 * D
		theta = 0.8 * Q;
	end	
	#correction de powel
	gammaK = (1-theta) * M * delta + theta * gammaL;

	if true ||(options.verb == 2) ##=== Impression ===##
		fprintf('--------------------------------BFGS------------------------------------------------------\n');
		fprintf('vap(M): ')
		fprintf('%f ',eig(M)')
		fprintf('\n%d \t theta = %f\n', info.niter, theta);
		fprintf('%12s %10s %10s %17s %15s %10s %10s \n',"g'd>0", "d'Md","gL'd", "d'Md/(d'Md-gL'd)", "gL'd/d'Md <0.2", "||gL||", "||d||");
		fprintf('%10.5e %10f %10f %17f %15f %10f %10f \n',gammaK'*delta,D,gammaL'*delta, D/D-S, S/D, norm(gammaL, inf), norm(delta, inf) );
	end ##============================## 
	
	eta = (gammaK'*gammaK) / (gammaK' * delta);
	if info.niter == 0 
		M = eta*eye(n);
	end
	
	M = M - (M*delta*delta' * M )/ (delta'*M*delta)  + eta;
	return
end