function [GR, WU, dd, ss, U, V, B, Q, QQ] = CMVToHessUnitaryReduction(A, U, V)
ss1=size(U);
n=ss1(1);
k=ss1(2);

S = zeros(2,0);
GR = zeros(2, 0);
WU  = zeros(3 * k, n);
dd = zeros(n, 1);
ss = zeros(n - 1, 1);

need_change_of_basis = (nargout == 9);

if need_change_of_basis
  QQ = eye(n);
end

for j = 1 : n - 2 * k - 1
    
    % Clean the j-th column
    for i = j + 2 * k - 1 : -1 : j + 2
        G = givens(A(i-1,j) + U(i-1,:) * V(j,:)', ...
                   A(i,j) + U(i,:) * V(j,:)');
        
        lo = max(1, i - 4*k - 2);
        hi = min(i + 4*k + 2, n);    
        
        A(i-1:i,lo:hi) = G * A(i-1:i,lo:hi);
        U(i-1:i,:) = G * U(i-1:i,:);
        
        S = PushRotation(S, G);
        
        if need_change_of_basis
          QQ(i-1:i,:) = G * QQ(i-1:i,:);
        end
        
        if j > 1
            GR = PushRotation(GR, G);
        end       
    end   
    
    % Store diagonal and subdiagonal entries
    dd(j) = A(j,j) + U(j,:) * V(j,:)'; 
    ss(j) = A(j+1,j) + U(j+1,:) * V(j,:)';
                
    % Now we can start applying the rotations that we have kept in standby...
    for i = j + 2 * k - 1 : -1 : j + 2
        [G, S] = PopRotation(S); 
                        
        A(:,i-1:i) = A(:,i-1:i) * G';
        V(i-1:i,:) = G * V(i-1:i,:);
        
        if (i > j + k)
            % Bulge chasing!
            s = i;
            
            while s <= n - 2 * k                
                G = givens(A(s-1+2*k,s-1), A(s+2*k,s-1));
                s = s + 2*k;
                
                lo = max(1, s - 4 * k - 1);
                hi = min(n, s + 4 * k + 1);
                
                A(s-1:s,lo:hi) = G * A(s-1:s,lo:hi);
                U(s-1:s,:) = G * U(s-1:s,:);
                V(s-1:s,:) = G * V(s-1:s,:);
                
                A(lo:hi,s-1:s) = A(lo:hi,s-1:s) * G';
                
                if need_change_of_basis
                  QQ(s-1:s,:) = G * QQ(s-1:s,:);
                end
            end
        end
                        
        if i == j + k + 1;
            S2 = zeros(2, 0);
            
            % Shift up the column before propagating it
            m = i;
            while m <= n
                for s = min(n, m + 2 * k - 1): -1 : m + k
                    G = givens(A(s-1,m-1), A(s,m-1));                    
                    
                    if m == i                        
                        lo = max(1, s - 3 * k - 1);
                        hi = min(n, s + 2 * k + 1);
                    else
                        lo = max(1, s - 2 * k - 1);
                        hi = min(n, s + 2 * k + 1);
                    end
                    
                    S2 = PushRotation(S2, G);
                    A(s-1:s,lo:hi) = G * A(s-1:s,lo:hi);
                    U(s-1:s,:) = G * U(s-1:s,:);
                    
                    if need_change_of_basis
                      QQ(s-1:s,:) = G * QQ(s-1:s,:);
                    end
                end
                m = m + 2 * k;
            end
            
            m = i;
            while m <= n
                for s = min(m + 2 * k - 1, n): -1 : m + k            
                    [G, S2] = PopRotation(S2);
                    
                    if m == i
                        lo = max(1, s - 3 * k - 1);
                        hi = min(n, s + 2 * k + 1);
                    else
                        lo = max(1, s - 2 * k - 1);
                        hi = min(n, s + 2 * k + 1);
                    end
                    
                    V(s-1:s,:) = G * V(s-1:s,:);                
                    A(lo:hi,s-1:s) = A(lo:hi,s-1:s) * G';
                    
                    if m == i && j > 1
                        GR = PushRotation(GR, G);
                    end
                end
                m = m + 2 * k;
            end
        end
    end
    
    % Here we store the diagonal and subdiagonal elements of the Hessenberg
    % form, and the j-th row from the superdiagonal element to position j +
    % 2 * k. This will be propagated in the next iterations by the
    % rotations which will be stored in GR. 
    l = min(3*k + 1, n - j);
    WU(1:l,j) = A(j,j+1:j+l)';
end


% Reduce the trailing matrix. 
B = A(n - 2*k : end, n - 2*k : end) + U(n-2*k:end,:) * V(n-2*k:end,:)';
[Q, B] = hess(B);

if need_change_of_basis
  QQ(n - 2*k : end, :) = Q' * QQ(n-2*k:end,:);
end

% In principle we would need to apply these transformations, but we are not
% doing since at the end we recover the bottom block from the Q and B. 
% 
% U(end-2*k:end,:) = Q' * U(end-2*k:end,:);
% V(end-2*k:end,:) = Q' * V(end-2*k:end,:);






