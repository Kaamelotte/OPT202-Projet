function [] = verifierGradient(xy)
  ##=== Variables utiles =====================================================##
  nn = length(xy)/2;
  x = xy(1:nn);
  y = xy(nn+1:end);
  [~,~,~,g,~,~,~,~] = chs(4,xy,[]);
  ##==========================================================================##
  
  ##=== Impression ===========================================================##
  fprintf('-----------------------------------------------\n');
  fprintf('%2s %6s %11s %9s %13s\n','i','pas',"f'(i)",'DF','erreur');
  fprintf('-----------------------------------------------\n');
  ##==========================================================================##
  
  DF = zeros(2*nn,1);
  err = 0;

  for i = 1:2*nn #pour chaque direction i
    e = zeros(2*nn,1); #vecteur de direction
    e(i) = 1;
    t = sqrt(eps) * max(1,xy(i)); #petit pas
    ##=== Calcul par différences finies ====================================##
    [ep,~,~,~,~,~] = chs(4, xy+t*e,[]);
    [em,~,~,~,~,~] = chs(4, xy-t*e,[]);
    DF(i) = ( ep - em )/(2*t);
    ##=== Comparaison à notre calcul du gradient ===========================##
    
    if DF(i) > 0
      err = abs(DF(i) - g(i) ) /DF(i);
      type = 'rel';
    else
      err = abs(DF(i) - g(i));
      type = 'abs';
    end
    fprintf('%2d %.2e %.5e %.5e %.5e %3s\n',i,t,g(i),DF(i),err,type);
  end
  fprintf('-----------------------------------------------\n');
  return
end