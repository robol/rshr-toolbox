function Y = MultiplyLeftSequence(Q, X)
%MULTIPLYLEFTSEQUENCE Y = Q * X

Y = X;
l = length(Q);
k = size(Q{1}) / 2;

for j = l : -1 : 1
    Y(end-(l-j+2)*k+1:end-(l-j)*k,:) = Q{j} * ...
        Y(end-(l-j+2)*k+1:end-(l-j)*k,:);
end


end

