% Test the backward stability of the reduction. 
% 
% In this example we test the backward error for different sizes and ranks.
% More precisely, we test the sizes
% 
%  N = [ 16, 32, 64, 128, 256, 512, 1024 ]
% 
% and the ranks 
% 
%  k = [ 1, 2, 4, 8, 16, 32 ]
%
% where some of them are tested only on "big enough" matrices. 

addpath ../matlab

Ns = [ 16, 32, 64, 128, 256, 512, 1024 ];
Ks = [ 2, 4, 8, 16, 32 ];

lk = length(Ks);
ln = length(Ns);

R = zeros(ln, lk);

for i = 1 : ln
    for j = 1 : lk
        n = Ns(i);
        k = Ks(j);
        
        if k > n / 4
            continue;
        end
        
        fprintf ('Testing n = %d, k = %d ...\n', n, k);
        
        [D,U,V] = rshr_build_example(n, k, 'ulr');
        [H, S] = rshr_ulr(D, U, V);
        H = rshr_full(H);
        R(i,j) = norm(H - S * (diag(D) + U*V') * S') / norm(diag(D + U*V'));
        
        fprintf ('N = %d, K = %d, BE = %e\n', R(i,j));
    end
end

dlmwrite('unitary_backward.dat', [ Ns', R ], '\t');
