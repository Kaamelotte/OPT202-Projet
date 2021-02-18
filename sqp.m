function [x,lm,info] = sqp(simul,x,lm,options)
################################################################################
## sqp :                                                                      ##
## Renvoie les multiplicateurs de Lagrange et les abscisses et ordonnées des  ##
## noeuds                                                                     ##
##                                                                            ##
## INPUT  - simul : spécification du simulateur                               ##
##        - x : vecteur contenant la valeur initiale                          ##
##        - lm : vecteur contenant les multiplcicateurs de Lagrange intiale   ##
##        - options : structure spéciafiant les paramètres de fonctionnement  ##
##                    de l’algorithme                                         ##
##                                                                            ##
## OUTPUT - x : vecteur contenant la valeur finale                            ##
##        - lm : vecteur contenant les multiplcicateurs de Lagrange finaux    ##                                                      ##
##        - info : structure donnant des informations sur le comportement du  ##
##                  solveur                                                   ##
##                                                                            ##
################################################################################
  
  ##=== Vérification de la consistance des arguments d'entrée ================##
  n = length(x);
  if mod(n,2) ~= 0
    info.status = 1 %length(x) impair
    return  
  end
  ##==========================================================================##
  
  ##=== Calcul des multiplicateurs lagrangiens ===============================##
  if length(lm) == 0
    [~,~,~,g,ae,ai,~,indic] = simul(4,x,lm);
    lm = -ae'\g;
  end
  m = length(lm);
  ##==========================================================================##
    
  ##=== Impression ===========================================================##
  if (options.verb == 1)
    fprintf('---------------------------------------------------------------------------------\n');
    fprintf('%4s %7s %10s %10s %10s %11s %9s %9s\n',...
            'iter','|gl|','|ce|','|x|','|lm|','alpha','phi','Q');
    fprintf('---------------------------------------------------------------------------------\n');
  end
  ##==========================================================================##

  ##=== Algorithme de Newton =================================================##
  nbIter = 0;
  nbSimul = 0;
  alpha = 1;  #pas initial = pas unité
  normFk = 1;
  while true
    
    ##=== Impression =========================================================##
    if options.verb == 1 
      fprintf('%4d ',nbIter);
    elseif options.verb == 2
      fprintf('-----------------------------------------------------\n');
      fprintf('iter %d, simul %d, ',nbIter,nbSimul);
    end
    ##========================================================================##
    
    [e,ce,ci,g,ae,ai,~,indic] = simul(4,x,lm);
    nbSimul += 1;
    grdl = g + ae' * lm;
    
    normFkp = max(norm(ce,Inf),norm(grdl,Inf));
    Q = normFkp / normFk^2;
    normFk = normFkp;
    
    ##=== Test d'optimalité ==================================================##
    if (norm(grdl,inf) < options.tol(1)) && (norm(ce,inf) < options.tol(2)) 
      info.status = 0; ##Solution trouvée
      info.niter = nbIter;
      ##=== Impression =======================================================##
      if options.verb == 1
        fprintf('%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',...
                norm(grdl,inf),norm(ce,inf),norm(x,inf),norm(lme,inf),alpha,phi,Q);
      end
      ##======================================================================##
      
      break;
    end
    ##========================================================================##
    
    ##=== Calcul de la direction de descente =================================##
    [~,~,~,~,~,~,hl,~] = simul(5,x,lm);
    nbSimul += 1;
    dF = [ hl, ae'; ae, zeros(m,m) ];
    F = [ g ; ce ];
    dir = -dF\F;
    
    ##========================================================================##
    
    phi = 0.5 * F' * F; #0.5*||F(z)||^2
    dphi = F' * dF;     #F(z)^T*F'(z)
    
    ##========================================================================##
    
    
    ##=== Impression =========================================================##
    if options.verb == 1
      fprintf('%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',...
              norm(grdl,inf),norm(ce,inf),norm(x,inf),norm(lme,inf),alpha,phi,Q);
    end
    ##========================================================================##
    
    ##=== Recherche linéaire pour le pas alpha ===============================##
    if options.rl == 0
      w = 1e-4;
      i = 0;       #nb d'itérations de la RL
      
      pente = w * dot(dphi,dir); #w.alpha_k.phi'(z_k).p_k
      
      ##=== Impression =================================================##
      if options.verb == 2
        fprintf('phi %.5e, pente %.5e\n\n',phi,pente);
        fprintf("Recherche linéaire d'Armijo: |d| = %.2e\n",norm(dir,inf))
        fprintf('%10s %15s %13s\n','alpha','phip-phi','DF(phi)');
      end
      ##================================================================##
      
      ##=== Boucle sur i =====================================================##
      while true
        xp = x + alpha*dir(1:n); 
        #lp = 
        lp = (1-alpha)*lme + alpha*dir(n+1:length(dir));

        
        [~,cep,cip,gp,aep,aip,~,~] = simul(4,xp,lp);
        nbSimul += 1;
        Fp = [gp; cep];
        phip = 0.5 * Fp' * Fp; #0.5*||F(z_k+alpha_k.p_k)||^2
        
        if options.verb == 2
          fprintf('%12.4e %14.5e %14.5e\n',alpha,phip-phi,(phip-phi)/alpha);
        end
        
        grdlp = gp + aep' * lmep;
        ##=== Test d'optimalité sur phi ======================================##
        if phip <= phi + alpha*pente
          if options.verb == 1
            fprintf('%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',...
                    norm(grdlp),norm(cep),norm(xp),norm(lmep),alpha,phi,Q);
          elseif options.verb == 2
            fprintf('|gl| = %.3e, |ce| = %.3e\n',norm(grdlp,inf),norm(cep,inf));
          end
          break;
        end
        ##====================================================================##
        
        alpha = alpha/2;
        i = i+1;
        ##===Test du nombre d'itérations déjà effectuées======================##
        if(i > 10)
            info.status = 3; #PB dans la RL
            if options.verb == 2
              fprintf('|gl| = %.3e, |ce| = %.3e\n',norm(grdlp,inf),norm(cep,inf));
            end
            break; 
        end
        ##====================================================================##
      end
      ##=== Fin boucle i =====================================================##
    end 
    ##=== Fin recherche linéaire =============================================##
    
    ##=== Calcul des nouveaux paramètres pour Newton =========================##
    x = x + alpha*dir(1:n);
    lme = (1-alpha)*lme + alpha*dir(n+1:length(dir));
    
    nbIter = nbIter + 1;
    info.niter = nbIter;
    ##========================================================================##
    
    ##===Test du nombre d'itérations déjà effectuées==========================##
    if(nbIter > options.maxit)
      info.status = 2; #Sortie de l'algo car pas de convergence
      break; 
    end
    ##========================================================================##
  end
  
  ##=== Impression ===========================================================##
  if options.verb == 1
    fprintf('---------------------------------------------------------------------------------\n');
  elseif options.verb == 2
    fprintf('-----------------------------------------------------\n');
  end
  ##==========================================================================##
  
  return
end
