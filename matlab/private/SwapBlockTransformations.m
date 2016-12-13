function [D, E, F] = SwapBlockTransformations(A, B, C)
%
% Given a sequence of block transformations of the form 
%
% [ A1 A2 0 ] [ I 0  0  ] [ C1 C2 0 ]
% [ A3 A4 0 ] [ 0 B1 B2 ] [ C3 C4 0 ] = Z
% [ 0  0  I ] [ 0 B3 B4 ] [ 0  0  I ]
%
% this function construct another sequence of transformations D, E, F such that
% 
% [ I 0  0  ] [ E1 E2 0 ] [ I 0  0  ]
% [ 0 D1 D2 ] [ E3 E4 0 ] [ 0 F1 F2 ] = Z
% [ 0 D3 D3 ] [ 0  0  I ] [ 0 F3 F4 ]
%
% All the above transformations are unitary. 

k = size(A, 2) / 2;

% Construct the matrix Z
Z = eye(3 * k);

Z(1:2*k,:) = C * Z(1:2*k,:);
Z(k+1:end,:) = B * Z(k+1:end,:);
Z(1:2*k,:)   = A * Z(1:2*k,:);

% Compute the matrices D, E and F
[D,~] = qr(Z(k+1:end,1:k));
Z(k+1:end,:) = D' * Z(k+1:end,:);

[E,~] = qr(Z(1:2*k,1:k));
Z(1:2*k,:) = E' * Z(1:2*k,:);

[F,~] = qr(Z(k+1:end,k+1:2*k));
Z(k+1:end,:) = F' * Z(k+1:end,:);

F = F * Z(k+1:end,k+1:end);
E(:,1:k) = E(:,1:k) * Z(1:k,1:k);