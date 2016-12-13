function Q = SequenceTimesDiagonal(Q, D)
%SEQUENCETIMESDIAGONAL Construct a new sequence by right multiplying by D. 

l = length(Q);
k = size(Q{1}) / 2; 

Q{l} = Q{l} * D(end-2*k+1:end,end-2*k+1:end);

for i = l-1 : -1 : 1
    Q{i}(:,1:k) = Q{i}(:,1:k) * D((i-1)*k+1:i*k, (i-1)*k+1:i*k);
end

end

