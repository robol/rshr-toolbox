function Y = MultiplyCMVSequence(Q, X)
%MULTIPLYCMVSEQUENCE Y = Q * X

Y = X;
l = length(Q);
k = size(Q{1}) / 2;

for j = l - 1 : -2 : 1
    Y((j-1)*k+1:(j+1)*k,:) = Q{j} * Y((j-1)*k+1:(j+1)*k,:);
end

for j = l : -2 : 1
    Y((j-1)*k+1:(j+1)*k,:) = Q{j} * Y((j-1)*k+1:(j+1)*k,:);
end

end

