function H = rshr_full(Hs)
%RHSH_FULL constructs the dense version of the Hessenberg matrix. 
% 
% H = RSHR_FULL(Hs) constructs the dense version of the structured
%       Hessenberg matrix encoded in the structure Hs. This structure
%       is expected to be as returned by one call to RSHR_DLR() or 
%       to RSHR_ULR(). 

structure_type = Hs.structure; 

if strcmp(structure_type, 'dlr')
    H = build_full_hess(Hs.dd, Hs.ss, Hs.U, Hs.V);
elseif strcmp(structure_type, 'ulr')
    H = BuildUpperHess(Hs.dd, Hs.ss, Hs.GR, Hs.WU, Hs.U, Hs.V, Hs.B, Hs.Q);
else
    error('Invalid structure type in the structure, aborting');
end

end
