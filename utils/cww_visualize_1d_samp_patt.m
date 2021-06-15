% Converts a one-dimensional sampling pattern to a matrix for visualization.
%
% This function takes a sampling pattern idx \subset {1,...,N}, and creates a
% matrix of size [round(0.15*N), N], where the columns indexed by `idx` is given 
% the value uint8(255), and columns not indexed by `idx` is given the value 0. 
% This matrix is ideal for visualizing a given sampling pattern.
%
% Arguments
% ---------
% N (int): Number of columns of the matrix and an upper bound for the elements 
%          in `idx`.
% idx (vector): Sampling pattern, given as a vector with integers in the range 
%               {1, ..., N}.
%
% Return
% ------
% X (mat, uint8): Matrix visualizing the sampling pattern `idx`.
function X = cww_visualize_1d_samp_patt(N, idx)

    X = zeros([round(0.15*N), N], 'uint8');
    X(:,idx) = uint8(255);

end
