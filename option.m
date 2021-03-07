function options = option(maxit = 50, quad = 1, verb = 1, deriv = 1, rl = 1, tol = [1.e-5,1.e-5,1.e-5], saveFig = 1, saveLog = 1 )
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
	options.saveFig = saveFig; #enregistrement sous jpg du res
	options.saveLog = saveLog; #enregistrement sous txt du res
	
	return
end;