function [SS, U, V, QQ] = Triangularize(SS, U, V)
%TRIANGULARIZE 

need_change_of_basis = (nargout == 4);

k = size(U, 2);
n = size(U, 1);

if need_change_of_basis
    QQ = eye(n);
end

SS(k+1:2*k,1:k) = triu(SS(k+1:2*k,1:k));

for j = 3 : 2 : n/k - 1
    [Q,R] = qr(SS(j*k+1:(j+1)*k,(j-2)*k+1:(j-1)*k));
    SS(j*k+1:(j+1)*k,:) = Q' * SS(j*k+1:(j+1)*k,:);
    SS(j*k+1:(j+1)*k,(j-2)*k+1:(j-1)*k) = R;
    
    SS(:,j*k+1:(j+1)*k) = SS(:,j*k+1:(j+1)*k) * Q;
    U(j*k+1:(j+1)*k,:) = Q' * U(j*k+1:(j+1)*k,:);
    V(j*k+1:(j+1)*k,:) = Q' * V(j*k+1:(j+1)*k,:);
    
    if need_change_of_basis
        QQ(j*k+1:(j+1)*k,:) = Q' * QQ(j*k+1:(j+1)*k,:);
    end
end

for j = 2 : 2 : n/k - 1
    [Q,R] = qr(SS((j-2)*k+1:(j-1)*k,j*k+1:(j+1)*k)');    
    SS(:,j*k+1:(j+1)*k) = SS(:,j*k+1:(j+1)*k) * Q;
    SS((j-2)*k+1:(j-1)*k,j*k+1:(j+1)*k) = R';
    SS(j*k+1:(j+1)*k, :) = Q' * SS(j*k+1:(j+1)*k, :);
    V(j*k+1:(j+1)*k,:) = Q' * V(j*k+1:(j+1)*k,:);
    U(j*k+1:(j+1)*k,:) = Q' * U(j*k+1:(j+1)*k,:);
    
    if need_change_of_basis
        QQ(j*k+1:(j+1)*k,:) = Q' * QQ(j*k+1:(j+1)*k,:);
    end
end


end

