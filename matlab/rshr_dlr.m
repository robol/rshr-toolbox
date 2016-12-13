function [H, P] = rshr_dlr(D, U, V)
%RSHR_DLR Hessenberg reduction for diagonal plus low rank matrices. 
%
% H = RSHR_DLR(D,U,V) performs a Hessenberg reduction of the
%        matrix A = diag(D) + U*V'. In particular, the input 
%        parameters are the following:
%
%        D: a vector of length n.
%        U: A n x k matrix.
%        V: A n x k matrix.
%
% [H, P] = RSHR_DLR(D, U, V) also returns the change of basis P such that 
%        P * (diag(D) + U*V') * P' is equal to H. Notice that in this case
%        the computation of P is only implemented in the MATLAB version, so
%        the MEX file will never be called. 
%
% The output is the Hessenberg matrix encoded in a structured form. The full
% Hessenberg matrix can be obtained by calling RSHR_FULL(H);

  % A = diag(D);
  n = length (D);
  k = size(U,2);

  n == size(U,1) || error ('U must have n rows where n = length(D)');
  n == size(V,1) || error ('V must have n rows where n = length(D)');
  k == size(V,2) || error ('V must have the same number of columns of U');

  dd = zeros (1, n);
  ss = zeros (1, n-1);

  % mex_available = exist ('hess_red2_impl', 'file') == 3;
  compute_change_of_basis = nargout == 2;

  if compute_change_of_basis
      %if ~mex_available
      %  warning ('MEX_NOT_FOUND', ...
      %  'hess_red2_impl MEX file not available, relying on the slow version');
      %end
      
      A = diag(D);

      P = eye(n);

      for i = n : -1 : 2
        for j = 1 : min (k, n - i + 1)
          l = i + j - 1;
          lo = max(1, l - k - 1);
          hi = min(l + k  +1, n);
          G = planerot (U(l-1:l,j));
          U(l-1:l,:) = G * U(l-1:l,:);
          V(l-1:l,:) = G * V(l-1:l,:);
          A(l-1:l,lo:hi) = G * A(l-1:l,lo:hi);
          A(lo:hi,l-1:l) = A(lo:hi,l-1:l) * G';
          P(l-1:l,:) = G * P(l-1:l,:);
        end

          % Perform bulge chasing if needed
          for s = i + k : n

            lo = max(1, s - k - 1);
            hi = min(s + k  +1, n);
            G = planerot (A(s-1:s,s-k-1));
            U(s-1:s,:) = G * U(s-1:s,:);
            V(s-1:s,:) = G * V(s-1:s,:);
            A(s-1:s,lo:hi) = G * A(s-1:s,lo:hi);
            A(lo:hi,s-1:s) = A(lo:hi,s-1:s) * G';
            P(s-1:s,:) = G * P(s-1:s,:);
         end
      end

      % Set zeros on the subdiagonal terms of U. This is done just to
      % make the output more readable.
      U = triu(U);

      % Step 2: Perform the actual Hessenberg reduction.
      for j = 1 : n - 2
        for i = min(n, j + k) : -1 : j + 2
          lo = max(1, i-k-1);
          hi = min(n, i+k+2);

          G = planerot (A(i-1:i,j) + U(i-1:i,:) * V(j,:)');
          U(i-1:i,:) = G * U(i-1:i,:);
          V(i-1:i,:) = G * V(i-1:i,:);
          A(i-1:i,lo:hi) = G * A(i-1:i,lo:hi);
          A(lo:hi,i-1:i) = A(lo:hi,i-1:i) * G';
          P(i-1:i,:) = G * P(i-1:i,:);

          for s = i + k : k : n
            lo = max (1, s - k - 1);
            hi = min (n, s+k+2);
            G = planerot (A(s-1:s,s-k-1) + U(s-1:s,:) * V(s-k-1,:)');
            U(s-1:s,:) = G * U(s-1:s,:);
            V(s-1:s,:) = G * V(s-1:s,:);
            A(s-1:s,lo:hi) = G * A(s-1:s,lo:hi);
            A(lo:hi,s-1:s) = A(lo:hi,s-1:s) * G';
            P(s-1:s,:) = G * P(s-1:s,:);
          end

          dd(j) = A(j,j) + U(j,:) * V(j,:)';
          ss(j) = A(j+1,j) + U(j+1,:) * V(j,:)';
        end 
      end

      dd(n-1:n) = diag(A(n-1:n,n-1:n) + U(n-1:n,:) * V(n-1:n,:)') ;
      ss(n-1) = A(n,n-1) + U(n,:) * V(n-1,:)';
  else
     % This version of hess_red2 is implemented in a relatively optimized
     % MEX file. It does not check the correctness of the arguments at all,
     % so they have to be checked in this function before passing them
     % to the MEX file.
     A = zeros(n,2*k+3);
     A(:,k+2) = D;
     [~,U,V,dd,ss] = hess_red2_wrapper (A, U, V);
   end

   H = struct (...
       'structure', 'dlr', ...
       'U', U, ...
       'V', V, ...
       'dd', dd, ...
       'ss', ss ...
   );
end
