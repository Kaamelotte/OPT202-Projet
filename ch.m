<<<<<<< Updated upstream
clear(); close(); clc();
global A B L R S color lgd;
color = 'r';  lgd = "resultat";
##=================== Cas test: ==============================================##
## TP1: '1 
## TP2:        '2a' '2b' '2c' '2d'
## TP3: '2d' '3a' '3b' '3c'
## TP4:        '4a' '4b' '4c'
## TP5:        '5a' '5b' '5c'
##=========================================================================##

castest = '2d';
test = '0';

##================ Optimiseur ===============================================##
options = option();
options.quad = 0;
options.verb = 2;
options.save = 0;
options.rl = 0;

if test == '0'	
	[x, lme, lmei, info] = res(castest, options, color_res_1 = [0.6,0.3,0.6]);
	info.status
=======
clear()
close()
clc()

global A B L R S color lgd

##=================== Cas test: ==============================================##
## TP1: '1																	  ##
## TP2: '2a' '2b' '2c' '2d'                                                   ##
## TP3: '2d' '3a' '3b' '3c'                                                   ##
## TP4: '4a' '4b' '4c'                                                    ##
##============================================================================##

castest = '4b'
color = 'r';
lgd = "resultat";
test = '0'

##================ Optimiseur ================================================##
options = option();

if test == '0'	
	[x, lme, lmi, info] = res(castest, options, color_res_1 = [0.6,0.3,0.6]);
	
>>>>>>> Stashed changes
	if info.status == 2
		fprintf('--------Condition non remplis--------\n');
		fprintf('%8s %10s %10s \n',...
					'|gl|','|ce|','|lmei-ci|');
		fprintf('%4e %10.4e %10.4e \n',...
					options.tol(1),options.tol(2), options.tol(3));
		fprintf('%4e %10.4e %10.4e \n',...
					info.tol(1),info.tol(2), info.tol(3));
		fprintf('\n--------------------------------------\n');
	end	
<<<<<<< Updated upstream
else##================ Simulateur ============================================##
	checker(castest, test);
end

=======
end


##================ Simulateur ================================================##
[L,xy,A,B,R,S] = casTest(castest);
switch test
	case '1' ##=== Trace de la chaine ===##	  
		disp('trace de la chaine');
		hold on;
		chs(1,xy,[], []);
		chs(6,xy,[], []);
		hold off;
	case 'c'  ##=== Energie potentielle et contraintes ======##
		disp("calcul de l'energie et des contraintes");
		[e,ce,ci,~,~,~,~,indic] = chs(2,xy,[],[])
	case 'g'  ##=== Gradient de e et jacobienne de c  =================##
		disp("calcul du gradient de e et de la jacobienne de c");
		[~,~,~,g,ae,ai,~,indic] = chs(4,xy,[],[])
	case 'hl' ##=== Hessien du lagrangien ============##
		disp("calcul du hessien du lagrangien");
		[~,~,~,~,~,~,hl,indic] = chs(5,xy,[],[]);
		hl = full(hl)
		valuerPropre = eig(hl)
		sum(eig(hl)>= 0)/length(eig(hl)) 
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
		E =  norm(hl - M,inf)
end##============================================================================##
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
>>>>>>> Stashed changes
if options.verb == 3 &&  test == '0'
	disp("Calcul de la hessienne reduite")
	[~,ce,ci,g,ae,ai,~,indic] = chs(4,x,lme,[]);
	[~,~,~,~,~,~,hl,indic] = chs(5,x,lme,[]);

	N = null(ae)
	eig(N'*hl*N)

	##=== Conditions suffisantes d'ordre 2 ===##
	fprintf('---------------------------------\n');
	fprintf('Verification des CS2:\n');
	fprintf('---------------------------------\n');
	fprintf('Grdl(x,lme) =\n');
	fprintf('%+25.5e\n',g+ae'*lme);
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
	    fprintf("|   <Grd^2l(x,lme).z|z> = %f\n",dot(hl*N(:,i),N(:,i)));
	    test += (dot(hl*N(:,i),N(:,i))<=0);
	end  
	fprintf('----------------------------------\n');
	
	fprintf('%.5f,\n%.5f,',x(1:length(x)/2),x(length(x)/2+1:length(x)))

  fprintf('\n')
end;
##=========================================================================##
