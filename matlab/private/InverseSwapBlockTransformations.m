function [D,E,F] = InverseSwapBlockTransformations(A, B, C)
% Inverse of the function SwapBlockTransformations(). 
%

A = flipud(fliplr(A));
B = flipud(fliplr(B));
C = flipud(fliplr(C));

[D, E, F] = SwapBlockTransformations(A, B, C);

D = flipud(fliplr(D));
E = flipud(fliplr(E));
F = flipud(fliplr(F));