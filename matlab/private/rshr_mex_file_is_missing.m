function missing = rshr_mex_file_is_missing()
%RSHR_MEX_FILE_IS_MISSING Check if hess_red2_impl.mex exists.
%
%
    path =  strcat(fileparts(which('rshr_dlr')), '/private/hess_red2_impl.mex');
    f = fopen(path, 'r');

    if f == -1
        missing = true;
    else
        missing = false;
        fclose(f);
    end
end
