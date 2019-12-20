% Two dimensional fast Walsh-Hadamard transform. The transform is normalized,
% so that the transform becomes unitary. 
%
% INPUT
% X     - Matrix to be transformed 
% order - The order of the Walsh-Hadamard transform. 
%           * 'sequency' (default)
%           * 'hadamard' 
%           * 'dyadic'
%
% OUTPUT
% Y - The transformed matrix
function Y = cww_fastwht_2d(X, order)

    if (nargin == 1)
        order = 'sequency';
    end
    N = size(X,1);
    if (exist('fastwht') == 3) % fastwht is installed
        Y = fastwht(X, [], order);
        Y = fastwht(Y', [], order)';
    else
        Y = fwht(X, [], order);
        Y = fwht(Y', [], order)';
    end
end


