function H = BuildUpperHess(dd, ss, GR, W, U, V, B, Q)
%BUILDUPPERHESS Construct the full upper Hessenberg form. 
%
% H = BuildUpperHess(DD, SS, GR, W) constructs the full upper Hessenberg 
% matrix defined by the parameters DD, SS, GR and W which have been 
% obtained by a call to {\ref UnitaryReduction2}. 
%
% Author: Leonardo Robol <leonardo.robol@cs.kuleuven.be>

% Recover the parameter of the rank and the dimensions from the input data.
n = length(dd);
k = (size(W, 1) - 1) / 3;

H = zeros(n);

for j = 1 : n - 2*k - 1
    l = min(3*k+1, n - j);
    H(j,j+1:j+l) = W(1:l,j)';
    
    if j < n -2*k - 1
        for s = j + 2 * k : -1 : j + 3
            [G, GR] = PopRotation(GR);
            H(:,s-1:s) = H(:,s-1:s) * G';        
        end
        
        for s = j + l : -1 : j + 2*k + 2
            [G, GR] = PopRotation(GR);
            H(:,s-1:s) = H(:,s-1:s) * G';
        end
    end
end

H = H + diag(dd) + diag(ss, -1) + triu(U * V', 1);
H(:, end - 2 * k : end) = H(:, end - 2 * k : end) * Q;
H(end-2*k : end,end-2*k : end) = B;

end

