function [x, lme, lmi, info, x2] = res(castest, options, ...
													  color_res_1 = [0.83,0.35,0.17], color_init = [0,0.42,0.7], color_res_2 = [0.6,0.3,0.6])

	global A B L R S color
	
	##=== Verification de la consistance des arguments d'entree ========================##
	if ( 0 > options.quad || options.quad > 3)
		info.status = 1 % non choix de l'algorithme
		return 
	end ##====================================================================##
	
	[L, xy, A, B, R, S] = casTest(castest);
	##=== Graphe de l'initialisation ===========================================##
		figure('Name',['Methode de Newton: cas test ',castest]);
		hold on;
		chs(6, xy );
		color = color_init;
		chs(1, xy );
		##=== Resolution par l'optimiseur ==========================================##
		if (options.quad == 0)
			[x,lme, lmi, info] = sqp(@chs, xy, [], [], options);
		else
			[x, lme, lmi, info] = oqs(@chs, xy, [], [], options);
		end;
		
		##=== Graphe de la solution ================================================##
		color = color_res_1;
		chs(1, x, lme, lmi);
		
		if(options.quad == 2)
			[x, lme, lmi, info] = oqs(@chs, xy, [], [], options);
			color = color_res_2;
			chs(1, x, lme, lmi);
			legend("", 'Initialisation','resultat_1', 'resultat_2');
		else
			legend('Initialisation','resultat');
		end
		
		title({ ['Methode de Newton: cas test  ',castest]; ...
			['Nbr Iterations: ', num2str(info.niter)]; "";"" });

		print(["figure_", castest,"_",num2str(options.rl), ".jpg"]);  
		hold off
		##==========================================================================##		
	L
  return
end