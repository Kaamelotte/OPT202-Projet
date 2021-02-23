function options = option(maxit = 10, quad = 1, verb = 1, deriv = 2, rl = 0, tol = [1.e-2,1.e-2,1.e-2], save = 1 )
	options.tol(1) = tol(1); # sur le gred du laplacien
	options.tol(2) = tol(2); # sur les conditions d egalite
	options.tol(3) = tol(3); # sur le min des multi de lagrange - les conditions d inegalite
	
	options.maxit = maxit;

	options.quad = quad;# Utilisation solveur quadratique ou non, obligatoire pour TP4
		#0: sans solveur quadratique
		#1: avec solveur quadratique
		#2:  avec et sans solveur quadratique

	options.deriv = deriv;
		#1: methode de quasi-Newton 
		#2: methode de Newton 

	options.rl = rl; #Recherche lineaire
	options.verb = verb;
	options.save = save; #enregistrement sous jpg du res
	
	return
end;