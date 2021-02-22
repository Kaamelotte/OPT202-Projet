function [e,ce,ci,g,ae,ai,hl,indic] = chs(indic,xy,lme,lmi)
################################################################################
## ch :
## Renvoie les multiplicateurs de Lagrange et les abscisses et ordonnees 
##                                                                   des noeuds
## INPUT  - indic : pilote du comportement
##        - xy : vecteur contenant les abscisses des ordonnees des noeuds
##        - lme : vecteur contenant les multiplcicateurs de Lagrange
## 
## OUTPUT - e : energie potnetielle
##        - ce : contraintes sur les longueurs des barres
##        - g : gradient de e
##        - ae : jacobienne des contraintes d'egalite
##        - hl : hessien du lagrangien
##        - indic : booléan du bon déroulement de la simulation
##
################################################################################
  global A B L R S color lgd
  ##=== Test de consistance des arguments d'entrée ===========================##
  if mod(length(xy),2) ~= 0
      indic = -1 % length(xy) impair
      return  
  end
  ##==========================================================================##
  
  ##=== Variables utiles =====================================================##
  nn = length(xy)/2; % nombre de noeuds intérieurs
  nb = nn+1; % nombre de barres
  p = length(R);
  
  xp = [xy(1:nn); A]; # { x_1, ..., x_nn }
  xm = [0; xy(1:nn) ]; # {x_0, ..., x_{nn-1} }
  
  yp = [xy(nn+1:2*nn); B]; # { y_1, ..., y_nn }
  ym = [0; xy(nn+1:2*nn)]; # {y_0, ..., y_{nn-1} }
  
  ##==========================================================================##

  ##=== Tracé de la chaîne======================================================##
	if indic == 1 
		 ##=== Tracé de la chaîne   === #
		  plt = plot([ 0; xp ], [ 0; yp ],'-o','linewidth',2, 'color',color, 'DisplayName', lgd);
		  Longueur = sqrt((xp-xm).^2 + (yp-ym).^2);
		  dx = 1/2*(xp-xm);
		  dy = 1/2*(yp-ym);
		  for i = 1:length(L)
			  text (xm(i)+dx(i), ym(i)+dy(i)+0.1, num2str(Longueur(i),4), ...
			  'color',color, 'fontweight', 'bold');
		  end;
		#axis off;
		indic = 0;
		return
	end ##====================================================================##	 
	  
  ##=== Tracé du plancher==========================================## 
	if indic == 6
	 ##=== Tracé du plancher === #
	  #calcule des limites de la figure
		if length(R) != 0
			d = 1;
		  xmin = min([0;xp]) -d;
		  xmax = max([xm;A]) +d;
		  ymin = min([0;yp])-d;
		  ymax = max([ym;A])+d ;
		 #		A
		 #			x
		 #         |  \
		 # 		|	 \
		 #			|	   \
		 #			|		 xB
		 #			|		|	
		 #		C	x___x D	  
		  x = [ xmin ; xmax ; xmax ; xmin]; # ( A, B , C , D )
		 shift = [ 0; 0; 2; 2 ];
		 
		 for i = 1:length(R)
			y = (R(i) + S(i)*x )- shift;
			y(3) = y(4);
			fill(x,y,  [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8], 'HandleVisibility','off')
		 end;
		  
		  #axis ([xmin+0.5, xmax-0.5, ymin, ymax], "square");
		  #axis off;
	  end;
      indic = 0;
      return
	end ##====================================================================##

  ##=== Calcul de e, ce et ci ================================================##
  if (indic == 2) || (indic == 4)
    x = xy(1:nn);  
    y = xy(nn+1:2*nn);
    
    ##===Energie potentielle===##
    e = 1/2*dot(L, yp + ym);
   
    ##===Contraintes=================##
    ce = (xp-xm).^2 + (yp-ym).^2 - L.^2;

    R_ = reshape(repmat(R,1,nn)',p*nn,1);#(R1 ... R1, R2 ... R2 ... Rp ... Rp)
    S_ = reshape(repmat(S,1,nn)',p*nn,1);
    ci = R_ + repmat(x,p,1).*S_ - repmat(y,p,1);
    
  ##==========================================================================##
  
  ##=== Calcul de g, ae et ai ================================================##
    if indic == 4  
      ##=== Gradient de l'énergie potentielle ===##
      g = [ zeros(nn,1); (L(1:nb-1)+L(2:nb))/2 ];
      
      ##=== Jacobienne des contraintes ============================##
		aex = spdiags([ [ (xm - xp)(2:nb) ;0], [xp-xm] ], -1:0, nb, nn);
		aey = spdiags([ [ (ym - yp)(2:nb) ;0], [yp-ym] ], -1:0, nb, nn);
		ae = sparse(2*[aex,aey]);
      
		aiy =full( repmat(-spdiags(ones(nn,1),0,nn,nn),p,1));
		aix= full(spdiags(repmat(S,nn,1), [0:-nn: -p*nn], p*nn , nn));

		ai = [ aix , aiy ];
      
		indic = 0;
		return
    end
      indic = 0;
      return
  end %if indic == 2 || indic == 4
  ##==========================================================================##
  
  ##=== Calcul de hl =========================================================##
  if indic == 5
	##=== Calcul des multiplicateurs lagrangiens d'égalité=============================##
    if length(lme) == 0 || length(lmi) == 0
		[~,~,~,g,ae,ai,~,indic] = chs(4,xy,lme, lmi);
		lme = -ae'\g;
		lmi = -ai'\g;
    end
	m = length(lme);
	##=======================================================================##
    
    diagcentre_e = [lme(1:nb-1) + lme(2:nb);lme(1:nb-1)+lme(2:nb)];
    diaginf_e = -[lme(2:nb-1); 0 ;lme(2:nb-1);0];
    diagsup_e = -[0;lme(2:nb-1); 0 ;lme(2:nb-1)];
    
    hl = 2*spdiags( [diaginf_e diagcentre_e diagsup_e], [-1 0 1 ], 2*nn, 2*nn );
    
    indic += 0;
    return
  end %if indic == 5
  ##==========================================================================##
  indic = -1 % valeur inconnue de indic
  return
end