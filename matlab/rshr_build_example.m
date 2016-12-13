function [D, U, V] = rshr_build_example(n, k, structure)
%RSHR_BUILD_EXAMPLE Build a random diagonal plus low rank matrix. 

D = randn(n, 1);

if exist('structure', 'var') && strcmp(structure, 'ulr')
    D = exp(2i * pi * D);
end

U = randn(n, k);
V = randn(n, k);

end

