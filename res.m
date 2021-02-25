function [x, lme, lmi, info, x2] = res(castest, options, ...
													  color_res_1 = [0.83,0.35,0.17], color_init = [0,0.42,0.7], color_res_2 = [0.6,0.3,0.6])

	global A B L R S color lgd
	
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
		lgd = "Initialisation";
		chs(1, xy);
		#legend(plt, "Initialisation");
		##=== Resolution par l'optimiseur ==========================================##
		lmi = [];
		if (options.quad == 0)
			[x,lme, info] = sqp(@chs, xy, [], options);
		else
			[x, lme, lmi, info] = oqs(@chs, xy, [], [], options);
		end;
		
		##=== Graphe de la solution ================================================##
		
		if(options.quad != 2)
			lgd = "Resultat";
			color = color_res_1;
			chs(1, x, lme, lmi);
		else 
			color = color_res_2;
			lgd = "Resultat_2";
			chs(1, x, lme, lmi);
			[x, lme, info] = sqp(@chs, xy, [], options);
			lgd = "Resultat_1";
			color = color_res_1;
			chs(1, x, lme, lmi);
		end
		legend;
		
		title({ ['Methode de Newton: cas test  ',castest]; ...
			['Nbr Iterations: ', num2str(info.niter)]; "";"" });
		
		if options.save == 1 
			print(["figure_", castest,"_",num2str(options.rl), ".jpg"]);  
		end
		hold off
		##==========================================================================##		
	L
  return
end