function options = option(maxit = 10, quad = 1, verb = 1, rl = 1, tol = [1.e-8,1.e-8,1.e-8], save = 1 )
	options.tol(1) = tol(1); # sur le gred du laplacien
	options.tol(2) = tol(2); # sur les conditions d egalite
	options.tol(3) = tol(3); # sur le min des multi de lagrange - les conditions d inegalite
	
	options.maxit = maxit;

	options.quad = quad;# Utilisation solveur quadratique ou non, obligatoire pour TP4
		#0 sans solveur quadratique
		# 1 avec solveur quadratique
		# 2  avec et sans solveur quadratique

	options.rl = rl; #Recherche lineaire
	options.verb = verb;
	options.save = save; #enregistrement sous jpg du res
	
	return
end;