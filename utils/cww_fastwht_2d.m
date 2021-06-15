% Two-dimensional fast Walsh-Hadamard transform. 
%
% Given an input matrix X, this function computes H*X*H', where H is a Hadamard 
% matrix.
%
% The transform is normalized so that the underlying linear operator is unitary. 
%
% Arguments
% ---------
% X (mat): Input matrix to be transformed 
% order (str):  The order of the Walsh-Hadamard transform. 
%               * 'sequency' (default)
%               * 'hadamard' 
%               * 'dyadic'
%
% Return
% ------
% Y (mat): The transformed matrix
%
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


