function [GR, WU, dd, ss, U2, V2, B, Q, S] = HessenbergReduction(D, U, V, debug_mode)
%HESSENBERGREDUCTION Perform the Hessenberg reduction of D + U*V'

if ~exist('debug_mode', 'var')
    debug_mode = false;
end

need_change_of_basis = (nargout == 9);

n = size(U, 1);

% 
% PART 0: In case the user wants to see the debug, perform some additional
% testing. 
% 
if debug_mode
    [~,~,conds] = condeig(diag(D) + U*V');
    [ce] = eig(diag(D) + U*V', 'nobalance');
    [ce,I] = sort(abs((ce)));
    conds = conds(I);
end

%
% PART 1: From diagonal to CMV
%
if debug_mode
    fprintf (' >> Reducing A to CMV plus spike\n');
    tic;
end

if need_change_of_basis
    [QL,U2,V2,~,S] = BlockUnitaryReduction(diag(D), U, V);
else
    [QL,U2,V2] = BlockUnitaryReduction(diag(D), U, V);
end

SS = MultiplyCMVSequence(QL, eye(n));
if need_change_of_basis
    [SS, U2, V2, QQ] = Triangularize(SS, U2, V2);
    S = QQ * S;
else
    [SS, U2, V2] = Triangularize(SS, U2, V2);
end

% U = U2; V = V2;

% norm(S * V - V2)
% norm(SS + U2 * V2' - S * (diag(D) + U*V') * S')

if debug_mode    
    ttocmv = toc;
    
    %
    % Checking the result by comparing the eigenvalues. 
    % 
    c3 = sort(abs(eig(SS + U2*V2', 'nobalance')));
    e3 = norm((c3 - ce) ./ conds) / norm(ce);
    fprintf ('      Reduction completed in %e seconds\n', ttocmv);
    fprintf ('      Backward error on the eigenvalues: %e\n\n', e3);
end

%
% PART 2: Reduction from CMV plus spike to upper Hessenberg.
%
if debug_mode
    fprintf(' >> Performing reduction from CMV plus spike to upper Hessenberg\n');
    tic;
end

if need_change_of_basis
  [GR, WU, dd, ss, U2, V2, B, Q, QQ] = CMVToHessUnitaryReduction(SS, U2, V2);
  S = QQ * S;
else
  [GR, WU, dd, ss, U2, V2, B, Q] = CMVToHessUnitaryReduction(SS, U2, V2);
end

if debug_mode
    ttohs = toc;
    
    HH = BuildUpperHess(dd, ss, GR, WU, U2, V2, B, Q);
    beh = norm(HH - S * (diag(D) + U*V') * S'); 

    %
    % Checking the results
    %
    c4 = sort(abs(eig(BuildUpperHess(dd, ss, GR, WU, U2, V2, B, Q), 'nobalance')));
    e4 = norm((c4 - ce) ./ conds) / norm(ce);
    fprintf ('      Reduction completed in %e seconds\n', ttohs);
    fprintf ('      Backward error on the eigenvalues: %e\n', e4);
    fprintf ('      Backward error of the Hessenberg form: %e\n', beh / norm(HH));
end

end

