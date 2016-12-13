function H = build_full_hess(D,S,U,V)
  H = diag(D) + diag(S,-1); 
  T = H - tril(U*V'); 
  H = H + triu(T', 1); 
  H = H + triu(U*V', 1);
end
