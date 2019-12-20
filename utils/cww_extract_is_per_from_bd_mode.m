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
