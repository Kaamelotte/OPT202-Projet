      nn = 8
	  p = 4
	  
	  S = 1:p
	 
      aiy =full( repmat(-spdiags(ones(nn,1),0,nn,nn),p,1))
	  aix= full(spdiags(repmat(S,nn,1), [0:-nn: -p*nn], p*nn , nn))
	  
      ai = [ aix , aiy ]