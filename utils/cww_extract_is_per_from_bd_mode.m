% Returns true of bd_mode is 'per'.
%
% Arguments
% ---------
% bd_mode (str): String specifying the wavelet boundary handling.
%
% Return
% ------
% is_per (int): Return 1 if bd_mode == 'per', 0 otherwise. 
%
function is_per = cww_extract_is_per_from_bd_mode(bd_mode)
    is_per = 0;
    if ischar(bd_mode)
        if strcmpi(bd_mode,'per')
            is_per = 1;
        elseif strcmpi(bd_mode,'bd')
            is_per = 0;
        else
            disp('Unknown bd_mode, using boundary wavelets');
        end
    else 
        disp('bd_mode: Must be a string')
    end
end
