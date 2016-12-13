function EigenvalueTest(n, k, structure)
%EIGENVALUETEST Check the computed eigenvalues

addpath ../matlab

[D,U,V] = rshr_build_example(n, k, structure);

if strcmp(structure, 'dlr')
    H = rshr_full(rshr_dlr(D, U, V));
else
    H = rshr_full(rshr_ulr(D, U, V));
end

sorted_eigs = @(A) sort(abs(eig(A)));

fprintf ('Residue on the eigenvalues: %e\n', ...
    norm(sorted_eigs(diag(D) + U*V') - sorted_eigs(H)) ./ ...
    norm(eig(diag(D) + U*V')));

end

