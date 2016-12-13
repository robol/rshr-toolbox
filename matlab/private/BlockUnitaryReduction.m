function [QL, U, V, QR, S] = BlockUnitaryReduction(D, U, V)
%BLOCKUNITARYREDUCTION Reduction to CMV plus spike. 

[QL,U] = BlockQR(U);

if nargout == 5
    S = MultiplyLeftSequence(QL, eye(size(U, 1)));
    need_change_of_basis = true;
else
    need_change_of_basis = false;
end

QR = QL; 
V = MultiplyLeftSequence(QR, V);
QL = SequenceTimesDiagonal(QL, D);

k = size(U, 2);

for i = length(QL) : -1 : 1
    % First we need to move the rotation downwards until it reaches the 
    % bottom of the sequence.
    Q = QR{i}';
    for j = i : length(QL) - 1
        [Q, QL{j}, QL{j+1}] = SwapBlockTransformations(QL{j}, QL{j+1}, Q);
        V(j*k+1:(j+2)*k,:) = Q' * V(j*k+1:(j+2)*k,:);
        
        if need_change_of_basis
            S(j*k+1:(j+2)*k,:) = Q' * S(j*k+1:(j+2)*k,:);
            % S(:,j*k+1:(j+2)*k) = S(:,j*k+1:(j+2)*k) * Q;
        end
    end
    
    % When this point has been reached the transformation is acting at the
    % bottom of the blocks, so we can fuse it with the left ones. 
    QL{end} = QL{end} * Q;
end

% Final step of block CMV reduction
for i = length(QL) : -2 : 2
    for j = length(QL) - 1 : -2 : i
        V((j-1)*k+1:(j+1)*k,:) = QL{j} * V((j-1)*k+1:(j+1)*k,:);
        
        if need_change_of_basis
            S((j-1)*k+1:(j+1)*k,:) = QL{j} * S((j-1)*k+1:(j+1)*k,:);
            % S(:,(j-1)*k+1:(j+1)*k) = S(:,(j-1)*k+1:(j+1)*k) * QL{j};
        end
    end
    
    for j = length(QL) : -2 : i
        V((j-1)*k+1:(j+1)*k,:) = QL{j} * V((j-1)*k+1:(j+1)*k,:);
        
        if need_change_of_basis
            S((j-1)*k+1:(j+1)*k,:) = QL{j} * S((j-1)*k+1:(j+1)*k,:);
            % S(:,(j-1)*k+1:(j+1)*k) = S(:,(j-1)*k+1:(j+1)*k) * QL{j};
        end
    end
end