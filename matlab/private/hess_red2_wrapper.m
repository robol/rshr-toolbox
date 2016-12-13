function [A,U,V,dd,ss] = hess_red2_wrapper(A, U, V)
%HESS_RED2_WRAPPER Call hess_red2_impl, and compile it if needed

if exist('hess_red2_impl', 'file') ~= 3
    % We need to ensure to be in the right directory
    folder = fileparts(which(mfilename));
    oldf = cd(folder);

    mex -largeArrayDims -c cdslib.f90 ...
        FFLAGS='-fPIC -fdefault-integer-8'
    mex -largeArrayDims -c hr_impl_full_to_kh.f90 ...
        FFLAGS='-fPIC -fdefault-integer-8'
    mex -largeArrayDims hess_red2_impl.c hr_impl_full_to_kh.o ...
        cdslib.o -lmwlapack -lmwblas -lgfortran
    cd(oldf);
end

[A,U,V,dd,ss] = hess_red2_impl(A, U, V);


end

