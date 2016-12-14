function [A,U,V,dd,ss] = hess_red2_wrapper(A, U, V)
%HESS_RED2_WRAPPER Call hess_red2_impl, and compile it if needed

% MATLAB / Octave detection
is_matlab = (exist ('OCTAVE_VERSION', 'builtin') == 0);

if (is_matlab && exist('hess_red2_impl', 'file') ~= 3) || ...
   (~is_matlab && rshr_mex_file_is_missing())

    % We need to ensure to be in the right directory
    folder = strcat(fileparts(which('rshr_dlr')), '/private');
    oldf = cd(folder);

    if is_matlab
        mex -largeArrayDims -c cdslib.f90 ...
            FFLAGS='-fPIC -fdefault-integer-8'
        mex -largeArrayDims -c hr_impl_full_to_kh.f90 ...
            FFLAGS='-fPIC -fdefault-integer-8'
        mex -largeArrayDims hess_red2_impl.c hr_impl_full_to_kh.o ...
            cdslib.o -lmwlapack -lmwblas -lgfortran
    else
        mex -c cdslib.f90
        mex -c hr_impl_full_to_kh.f90
        mex hess_red2_impl.c hr_impl_full_to_kh.o cdslib.o ...
            -llapack -lblas -lgfortran

       % Workaround needed for Octave, we need to clear functions
       % from memory to allow finding the new one.
       clear -f
    end
    cd(oldf);
end

[A,U,V,dd,ss] = hess_red2_impl(A, U, V);


end

