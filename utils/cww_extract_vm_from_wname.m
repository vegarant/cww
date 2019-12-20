function vm = cww_extract_vm_from_wname(wname)
    
    vm = 4;
    if ischar(wname) && length(wname) > 2
        if strcmpi(wname(1:2), 'db')
            vm = str2num(wname(3:end));
        elseif strcmpi(wname(1:3), 'sym')
            vm = str2num(wname(4:end));
        elseif strcmpi(wname, 'haar')
            vm = 1;
        else
            sprintf('wname must be on the form dbX or symX where X is an integer,\n')
            sprintf('vm defaults to %d\n', vm);
        end
    else
        sprintf('wname must be a string, vm defaults to %d\n', vm);
    end
    if isempty(vm)
        sprintf('cww_extract_vm_from_wname failed to extract vm.\n');
        vm = 4;
        sprintf('vm defaults to %d\n', vm);
    end
end
