clear()
close()

global A B L R S color

##=================== Cas test: ==============================================##
## TP1: '1																	  ##
## TP2: '2a' '2b' '2c' '2d'                                                   ##
## TP3: '2d' '3a' '3b' '3c'                                                   ##
##============================================================================##

castest = '4b'
[L,xy,A,B,R,S] = casTest(castest);
color = 'r';
test = '0'

##================ Simulateur ================================================##
switch test
	case '1' ##=== Trace de la chaine ===##	  
	  disp('trace de la chaine');
	  hold on 
	  chs(1,xy,[], []);
	  chs(6,xy,[], []);
	  hold off
	case '2'  ##=== Energie potentielle et contraintes ======##
	  disp("calcul de l'energie et des contraintes");
	  [e,ce,ci,~,~,~,~,indic] = chs(2,xy,[])
	case '4'  ##=== Gradient de e et jacobienne de c  =================##
		disp("calcul du gradient de e et de la jacobienne de c");
		[~,~,~,g,ae,ai,~,indic] = chs(4,xy,[])
	case '5' ##=== Hessien du lagrangien ============##
		disp("calcul du hessien du lagrangien");
		[~,~,~,~,~,~,hl,indic] = chs(5,xy,[])
		full(hl)
	case 'grad' ##=== Verification du gradient de l'energie potentielle ===##
		disp("verification du gradient de l'energie potentielle"); 
		verifierGradient(xy);
end 

##============================================================================##

##================ Optimiseur ================================================##
options.tol(1) = 1.e-8; # sur le gred du laplacien
options.tol(2) = 1.e-8; # sur les conditions d egalite
options.tol(3) = 1.e-8; # sur le min des multi de lagrange - les conditions d inegalite
options.maxit = 10;

options.rl = 1; #Recherche lineaire
options.verb = 2;
if options.rl == 1 && options.verb == 2
	options.verb = 1;
end;

if test == '0'    
	##=== Graphe de l'initialisation ===========================================##
	figure('Name',['Methode de Newton: cas test ',castest]);
	color = [0,0.42,0.7];
	hold on;
	chs(1,xy,[]);
  
	##=== Resolution par l'optimiseur ==========================================##
	[x,lm,info] = sqp(@chs, xy, [], [], options);
    
	##=== Graphe de la solution ================================================##
	color = [0.83,0.35,0.17];
	chs(1,x,lm, []);
	chs(6,x,lm, []);
	legend('Initialisation','resultat');

	title({ ['Methode de Newton: cas test  ',castest]; ...
		['Nbr Iterations: ', num2str(info.niter)]; "";"" });

	print(["figure_", castest,"_",num2str(options.rl), ".jpg"]);  
	hold off
	##==========================================================================##
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	if options.verb == 3
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
	end;
	
	
  fprintf('%.5f, ',x(1:length(x)/2))
  fprintf('\n') 
  fprintf('%.5f, ',x(length(x)/2+1:length(x)))

  fprintf('\n')
end
##============================================================================##
