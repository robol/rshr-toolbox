% Perform some timing measurements. 

addpath ../matlab

k = 4;
n = 16;

nk = 9;
tk = zeros(1, nk);

for i = 1 : nk
    kk = 2^(i-1);
    nn = 8192;

    [D, U, V] = rshr_build_example(nn, kk, 'dlr');
    tk(i) = timeit (@()  rshr_dlr(D, U, V), 1);
    
    fprintf ('  Time for n = %d, k = %d, %e\n', nn, kk, tk(i));
end

dlmwrite('hermitian_tk.dat', [ 2.^(0:nk-1)' , tk' ], '\t');

% Test complexity in the size, given a fixed block size
nn = 12; 
tn = zeros(1, nn);
tnf = zeros(1, nn);

for i = 1 : nn
    k = 4;
    n = 2^(i+2);
    
    [D, U, V] = rshr_build_example(n, k, 'dlr');
    tn(i) = timeit (@()  rshr_dlr(D, U, V), 1);
    
    fprintf ('  Time for n = %d, k = %d, %e\n', n, k, tn(i));
    
    % We also compare with a full Hessenberg reduction
    tnf(i) = timeit (@() hess(diag(D) + U*V'), 1);
end

dlmwrite('hermitian_tn.dat', [ 2.^(3:nn+2)' , tn' ], '\t');
dlmwrite('hermitian_tnf.dat', [ 2.^(3:nn+2)' , tnf' ], '\t');