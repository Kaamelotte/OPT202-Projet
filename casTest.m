function [L,xy,A,B, R, S] = casTest(numcas)
	R = [];
	S = [];
  ##==================================TP1===================================##
  if numcas == '1' || numcas == '2a'
    #Longueurs des barres
    L = [0.7, 0.5, 0.3, 0.2, 0.5]';
    #Position initiale des noeuds
    xy = [0.2, 0.4, 0.6, 0.8,...
         -1.0, -1.5, -1.5, -1.3]';
    A = 1.0;  B = -1.0;

  ##==================================TP2===================================##
elseif numcas == '2b'
    #Longueurs des barres
    L = [0.7, 0.5, 0.3, 0.2, 0.5]';
    #Position initiale des noeuds 
    xy = [0.2, 0.4, 0.6, 0.8,...
          1.0, 1.5, 1.5, 1.3]';
    A = 1.0;  B = -1.0;
  elseif numcas == '2c'
    #Longueurs des barres
    L = [0.7, 0.5, 0.3, 0.2, 0.5]';
    #Position initiale des noeuds 
    xy = [0.2, 0.4, 0.6, 0.8,...
         -1.0, -1.5, 1.5, -1.3]';
    A = 1.0;  B = -1.0;
  elseif numcas == '2d'
    #Longueurs des barres
    L = [0.7, 0.5, 0.3, 0.2, 0.5]';
    #Position initiale des noeuds 
    xy = [0.2, 0.4, 0.6, 0.8,...
         1.0, -1.2, 1.5, -1.3]';
    A = 1.0;  B = -1.0;
  
  ##==================================TP3===================================##
  elseif numcas == '3a'
    #Longueurs des barres
    L = [0.6, 0.6]';
    #Position initiale des noeuds 
    xy = [0.5, ...
          0.4]';
    A = 1.0;  B = 0.0;
  elseif numcas == '3b'
    #Longueurs des barres
    L = [2.0, 1.0]';
    #Position initiale des noeuds 
    xy = [0.5,...
          0.3]';
    A = 1.0;  B = 0.0;
  elseif numcas == '3c'
    #Longueurs des barres
    L = [2.0, 1.0]';
    #Position initiale des noeuds 
    xy = [0.3,...
          0.3]';
    A = 0.0;  B = -1.0;
  ##==================================TP4===================================##
  elseif numcas == '4a'
    #Longueurs des barres
    L = [0.7 0.5 0.3 0.2  0.5]';
    #Position initiale des noeuds 
    xy = [ 	0.2 0.4 0.6 0.8 ...
				1 1.5 1.5 1.3]';

    A = 1.0;  B = 0.0;
  elseif numcas == '4b'
    #Longueurs des barres
    L = [0.2 0.2 0.2 0.3 0.3 0.5 0.2 0.2 0.3 0.1]';
    #Position initiale des noeuds 
    xy = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
			-0.5 -0.9 -1.2 -1.4 -1.5 -1.4 -1.2 -0.9 -0.5]';

    A = 1.0;  B = 0.0;
	R = [-0.25]
	S = [-0.5]
  elseif numcas == '4c'
    #Longueurs des barres
    L = [0.2 0.2 0.2 0.3 0.3 0.5 0.2 0.2 0.3 0.1]';
    #Position initiale des noeuds 
    xy = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
			-0.5 -0.9 -1.2 -1.4 -1.5 -1.4 -1.2 -0.9 -0.5]';

    A = 1.0;  B = 0.0;
	R = [-0.25; -0.5];
	S = [-0.5; 0];

  end
  return
end