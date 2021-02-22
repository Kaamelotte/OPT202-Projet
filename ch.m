clear()
close()

global A B L R S color

##=================== Cas test: ==============================================##
## TP1: '1																	  ##
## TP2: '2a' '2b' '2c' '2d'                                                   ##
## TP3: '2d' '3a' '3b' '3c'                                                   ##
## TP4: '4a' '4b' '4c'                                                    ##
##============================================================================##

castest = '2a'
color = 'r';
test = '0'

##================ Optimiseur ================================================##
options.tol(1) = 1.e-8; # sur le gred du laplacien
options.tol(2) = 1.e-8; # sur les conditions d egalite
options.tol(3) = 1.e-8; # sur le min des multi de lagrange - les conditions d inegalite
options.maxit = 10;
options.quad = 2
	#0 sans solveur quadratique
	# 1 avec solveur quadratique
	# 2 avec et sans solveur quadratique

options.rl = 1; #Recherche lineaire
options.verb = 0; #Choix affichage

if test == '0'
	
	options.quad = 1
	[x, lme, lmi, ~] = res(castest, options);
	
	options.quad = 2
	[x, lme, lmi, ~] = res(castest, options);

end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	if options.verb == 3 &&  test == '0'
		##=== calcul de la hessienne reduite =======================================##
		disp("Calcul de la hessienne reduite")
		[~,ce,ci,g,ae,ai,~,indic] = chs(4,x,lm,[]);
		[~,~,~,~,~,~,hl,indic] = chs(5,x,lm,[]);

		N = null(ae)
		eig(N'*hl*N)

		##=== Conditions suffisantes d'ordre 2 =====================================##
		fprintf('---------------------------------\n');
		fprintf('Verification des CS2:\n');
		fprintf('---------------------------------\n');
		fprintf('Grdl(x,lm) =\n');
		fprintf('%+25.5e\n',g+ae'*lm);
		fprintf('---------------------------------\n');
		fprintf('c(x) =\n');
		fprintf('%+25.5e\n',ce);
		fprintf('---------------------------------\n');
		fprintf("hl defini positif sur\n le noyau de c'(x)?\n");
		fprintf('---------------------------------\n');
		%fprintf("Noyau de c'(x):\n");
		N = null(ae)
		%fprintf('---------------------------------\n');
		test = 0;
		for i=1:size(N,2)
		fprintf("|   <Grd^2l(x,lm).z|z> = %f\n",dot(hl*N(:,i),N(:,i)));
		test += (dot(hl*N(:,i),N(:,i))<=0);
		end  
		fprintf('----------------------------------\n');
		
		fprintf('%.5f,\n%.5f,',x(1:length(x)/2),x(length(x)/2+1:length(x)))

	  fprintf('\n')
	end;
end
##============================================================================##



##================ Simulateur ================================================##
[L,xy,A,B,R,S] = casTest(castest);
switch test
	case '1' ##=== Trace de la chaine ===##	  
		disp('trace de la chaine');
		hold on 
		chs(1,xy,[], []);
		chs(6,xy,[], []);
		hold off
	case '2'  ##=== Energie potentielle et contraintes ======##
		disp("calcul de l'energie et des contraintes");
		[e,ce,ci,~,~,~,~,indic] = chs(2,xy,[],[])
	case '4'  ##=== Gradient de e et jacobienne de c  =================##
		disp("calcul du gradient de e et de la jacobienne de c");
		[~,~,~,g,ae,ai,~,indic] = chs(4,xy,[],[])
	case '5' ##=== Hessien du lagrangien ============##
		disp("calcul du hessien du lagrangien");
		[~,~,~,~,~,~,hl,indic] = chs(5,xy,[],[])
		full(hl)
	case 'grad' ##=== Verification du gradient de l'energie potentielle ===##
		disp("verification du gradient de l'energie potentielle"); 
		verifierGradient(xy);
	case 'cholmod' ##=== Test de la fonction cholmod ===##
		disp("Test de la fonction cholmod\n"); 
		#Tolerance
		tol_small = 1.e-5; 
		tol_big = 1.e+5;
		[~,~,~,~,~,~,hl,indic] = chs(5,xy,[],[]);
		[L, d, flag] = cholmod(hl, tol_small, tol_big);
		full(hl)
		full(L)
		D = diag(d)
		M = L*D*L'
		disp("M est bien defini positive\n")
		eig(M) #valeur propre de M
		disp("Erreur d'approximation de hl\n")
		E =  sum(sum(abs(hl - M)) )
end##============================================================================##
