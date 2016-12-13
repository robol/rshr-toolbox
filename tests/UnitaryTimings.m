2% Perform some timing measurements. 

addpath ../matlab

k = 4;
n = 16;

nk = 5;
tk = zeros(1, nk);

for i = 1 : nk
    kk = 2^(i-1);
    nn = 8 * 2^nk;

    [D, U, V] = rshr_build_example(nn, kk, 'ulr');
    tk(i) = timeit (@()  rshr_ulr(D, U, V), 1);
    
    fprintf ('  Time for n = %d, k = %d, %e\n', nn, kk, tk(i));
end

dlmwrite('unitary_tk.dat', [ 2.^(0:nk-1)' , tk' ], '\t');

% Test complexity in the size, given a fixed block size
nn = 7; 
tn = zeros(1, nn);

for i = 1 : nn
    k = 4;
    n = 2^(i+3);
    
    [D, U, V] = rshr_build_example(n, k, 'ulr');
    tn(i) = timeit (@()  rshr_ulr(D, U, V), 1);
    
    fprintf ('  Time for n = %d, k = %d, %e\n', n, k, tn(i));
end

dlmwrite('unitary_tn.dat', [ 2.^(4:nn+3)' , tn' ], '\t');