function [x, lme, lmi, info, x2] = res(castest, options,color_res_1 = [0.83,0.35,0.17], color_init = [0,0.42,0.7], color_res_2 = [0.6,0.3,0.6])
    
	global A B L R S color lgd
	
	##=== Verification de la consistance des arguments d'entree ========================##
	if ( 0 > options.quad || options.quad > 3)
		info.status = 1 % non choix de l'algorithme
		return 
	end ##====================================================================##

	if options.saveLog == 1
		diary(["res_", castest,"_",num2str(options.rl), ".txt"]);
	end
	
	[L, xy, A, B, R, S] = casTest(castest);
	##=== Graphe de l'initialisation ===========================================##
    figure('Name',['Methode de Newton: cas test ',castest]);
    hold on;
    chs(6, xy );#dessine le plancher
    color = color_init;
    lgd = "Initialisation";
    chs(1, xy);

		##=== Resolution par l'optimiseur ==========================================##
		lmi = [];
	    [x,lme, lmi info] = sqp(@chs, xy, [], [], options);
        
        if options.verb == 5 #=== Impression ===#
	        disp("Calcul de la hessienne reduite")
	        [e,ce,ci,g,ae,ai,~,indic] = chs(4,x,lme,[]);
	        [~,~,~,~,~,~,hl,indic] = chs(5,x,lme,[]);
	        N = null([ae;ai])
	        vap = eig(N'*hl*N)
            fprintf('Verification des CS2:\n');
            fprintf('Grdl(x,lm) =\n');   fprintf('%+25.5e\n',g + [ae ; ai ]' * [lme ; lmi ]);     fprintf('\n');
            fprintf('ce(x) =\n');           fprintf('%+25.5e\n',ce);                                           fprintf('\n');
            fprintf('ci(x) =\n');            fprintf('%+25.5e\n',ci);                                            fprintf('\n');
            fprintf("hl defini positif sur\n le noyau de c'(x)?\n");
            N = null([ae;ai]);
            test = 0;
            for i=1:size(N,2)
                fprintf("|   <Grd^2l(x,lm).z|z> = %f\n",dot(hl*N(:,i),N(:,i)));
                test += (dot(hl*N(:,i),N(:,i))<=0);
            end  
            fprintf('----------------------------------\n');
            fprintf('%.5f,',x(1:length(x)/2));                    fprintf("\n");
            fprintf('%.5f,',x(length(x)/2+1:length(x)));    fprintf('\n');
        end
        if info.status == 2
            fprintf('--------Condition non remplie--------\n');
		    fprintf('%8s %10s %13s \n','|gl|','|ce|','|lmi-ci|');
		    fprintf('%4e %10.4e %10.4e \n',options.tol(1),options.tol(2), options.tol(3));
		    fprintf('%4e %10.4e %10.4e \n',info.tol(1),info.tol(2), info.tol(3));
		    fprintf('\n--------------------------------------\n');
        end
		##=== Graphe de la solution ================================================##
		figure(1);
		if (options.quad != 2)
			lgd = "Resultat";
			color = color_res_1;
			chs(1, x, lme, lmi);
		else 
			color = color_res_2;
			lgd = "Resultat_2";
			chs(1, x, lme, lmi);
			
			[x, lme, lmi, info] = sqp(@chs, xy, [], [], options);
			lgd = "Resultat_1";
			color = color_res_1;
			chs(1, x, lme, lmi);
		end
		legend;
		title({ ['Methode de Newton: cas test  ',castest];['Nbr Iterations: ', num2str(info.niter)]; "";"" });
		
		if options.saveFig == 1 
			print(["figure_", castest,"_",num2str(options.rl), ".jpg"]);  
		end
        diary off;
	##==========================================================================##		
	L
  return
end