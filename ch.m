clear(); close(); clc();
global A B L R S color lgd;
color = 'r';  lgd = "resultat";
##=================== Cas test: ==============================================##
## TP1: '1																	  ##
## TP2: '2a' '2b' '2c' '2d'                                                   ##
## TP3: '2d' '3a' '3b' '3c'                                                   ##
## TP4: '4a' '4b' '4c'                                                    ##
##============================================================================##

castest = '4b'
test = 'g'

##================ Optimiseur ================================================##
options = option();
options.save = 0;

if test == '0'	
	[x, lme, lmei, info] = res(castest, options, color_res_1 = [0.6,0.3,0.6]);
	
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
else##================ Simulateur ================================================##
	checker(castest, test);
end

if options.verb == 3 &&  test == '0'
	##=== calcul de la hessienne reduite =======================================##
	disp("Calcul de la hessienne reduite")
	[~,ce,ci,g,ae,ai,~,indic] = chs(4,x,lme,[]);
	[~,~,~,~,~,~,hl,indic] = chs(5,x,lme,[]);

	N = null(ae)
	eig(N'*hl*N)

	##=== Conditions suffisantes d'ordre 2 =====================================##
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
##============================================================================##
