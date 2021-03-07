clear(); close(); clc(); fclose('all');

global A B L R S color lgd;
color = 'r';  lgd = "resultat";
##=================== Cas test: ==============================================##
## TP1: '1 
## TP2:        '2a' '2b' '2c' '2d'
## TP3: '2d' '3a' '3b' '3c'
## TP4:        '4a' '4b' '4c'
## TP5:        '5a' '5b' '5c'
##=========================================================================##

castest = '5c';
test = '0';
options = option(100);
options.verb ==1

if test == '0'  ##=== Optimiseur ===##	
	[x, lme, lmi, info] = res(castest, options, color_res_1 = [0.6,0.3,0.6]);
else                ##=== Simulateur ===##
	checker(castest, test);
end