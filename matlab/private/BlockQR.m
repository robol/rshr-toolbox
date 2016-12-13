function [Q, U] = BlockQR(U)
%BLOCKQR Compute a sequence of transformations stored in the cell array Q
%that takes the matrix U in upper triangular form. 

k = size(U, 2);
n = size(U, 1);

Q = {};

for i = n - 2*k : -k : 0
    [QQ,R] = qr(U(i+1:i+2*k, :));
    Q = { QQ', Q{:} };
    U(i+1:i+2*k,:) = R;
end


end

