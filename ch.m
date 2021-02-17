clear()
close()
global A B L R S color

##=================== Cas test: ==============================================##
## TP1: '1'       test                                                            ##
## TP2: '2a' '2b' '2c' '2d'                                                   ##
## TP3: '2d' '3a' '3b' '3c'                                                   ##
##============================================================================##

castest = '2a'
[L,xy,A,B] = casTest(castest);
color = 'r';
test = 'grad'

##================ Simulateur ================================================##
if test == '1'
  ##=== Tracé de la chaine ===##
  disp('tracé de la chaine');
  chs(1,xy,[]);
elseif test == '2'
  ##=== Energie potentielle et contraintes ======##
  disp("calcul de l'énergie et des contraintes");
  [e,ce,ci,~,~,~,~,indic] = chs(2,xy,[])
elseif test == '4'
  ##=== Gradient de e et jacobienne de c  =================##
  disp("calcul du gradient de e et de la jacobienne de c");
  [~,~,~,g,ae,ai,~,indic] = chs(4,xy,[])
elseif test == '5'
  ##=== Hessien du lagrangien ============##
  disp("calcul du hessien du lagrangien");
  [~,~,~,~,~,~,hl,indic] = chs(5,xy,[])
  full(hl)
elseif test == 'grad'
  ##=== Vérification du gradient de l'énergie potentielle ===##
  disp("vérification du gradient de l'énergie potentielle"); 
  verifierGradient(xy);
  ##=========================================================##
end
##============================================================================##

##================ Optimiseur ================================================##
options.tol(1) = 1.e-8;
options.tol(2) = 1.e-8;
options.maxit = 100;
options.rl = 1;
options.verb = 1;
if test == '0'
    
  ##=== Graphe de l'initialisation ===========================================##
  figure('Name',['Méthode de Newton: cas test ',castest]);
  hold on;
  color = [0,0.42,0.7];
  chs(1,xy,[]);
  
  ##=== Résolution par l'optimiseur ==========================================##
  [x,lm,info] = sqp(@chs,xy,[],options);
    
  ##=== calcul de la hessienne réduite =======================================##
  disp("Calcul de la hessienne réduite")
  [~,ce,ci,g,ae,ai,~,indic] = chs(4,x,lm);
  [~,~,~,~,~,~,hl,indic] = chs(5,x,lm);
  
  N = null(a)
  eig(N'*hl*N)
    
  ##=== Conditions suffisantes d'ordre 2 =====================================##
  fprintf('---------------------------------\n');
  fprintf('Vérification des CS2:\n');
  fprintf('---------------------------------\n');
  fprintf('Grdl(x,lm) =\n');
  fprintf('%+25.5e\n',g+ae'*lm);
  fprintf('---------------------------------\n');
  fprintf('c(x) =\n');
  fprintf('%+25.5e\n',ce);
  fprintf('---------------------------------\n');
  fprintf("hl défini positif sur\n le noyau de c'(x)?\n");
  fprintf('---------------------------------\n');
  %fprintf("Noyau de c'(x):\n");
  N = null(a)
  %fprintf('---------------------------------\n');
  test = 0;
  for i=1:size(N,2)
    fprintf("|   <Grd^2l(x,lm).z|z> = %f\n",dot(hl*N(:,i),N(:,i)));
    test += (dot(hl*N(:,i),N(:,i))<=0);
  end  
  fprintf('----------------------------------\n');
  
  ##=== Graphe de la solution ================================================##
  color = [0.83,0.35,0.17];
  chs(1,x,lm);
  legend('Initialisation','résultat');
  if test == 0
    title({ ['Méthode de Newton: cas test  ',castest]; ...
            ['Nbr Iterations: ', num2str(info.niter)]; ...
            "hl définie positive sur Ker(c'(x))"; "";"" });
   
  else
    title({ ['Méthode de Newton: cas test  ',castest]; ...
            ['Nbr Iterations: ', num2str(info.niter)]; ...
            "hl non définie positive sur Ker(c'(x))"; "";"" });
  end
  
  print(["figure_", castest,"_",num2str(options.rl), ".jpg"]);  
  ##==========================================================================##
  fprintf('%.5f, ',x)
  fprintf('\n')
end
##============================================================================##
