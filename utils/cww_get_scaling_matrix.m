% Creates a matrix where each column contain the coefficients of exactly one
% scaling function. All, except the `vm` first and last columns, contains
% translated versions of the same scaling function. The returned matrix have
% size [N, M]. Thus this function is ideal if you need to plot a wavelet
% approximation to a function, given some wavelet coefficents
%
% INPUT
% vm - Number of vanishing moments for a Daubechies wavelet
% N - Number of row in output the matrix
% M - Number of columns in the output matrix
% It is assumed that M = 2^r and N = 2^(r+q) for positive integers r and q
%
% RETURNS
% A N × M matrix
function A = cww_get_scaling_matrix(log2N, log2M, wname, bd_mode)
    if (log2M + 1 > log2N)
        disp('Error: N should satisfy N >= 2*M');
    end
    is_per = cww_extract_is_per_from_bd_mode(bd_mode)
    if is_per
        A = cww_get_scaling_matrix_per(log2N, log2M, wname);
    else
        A = cww_get_scaling_matrix_bd(log2N, log2M, wname);
    end
end

% Creates a matrix where each column contain the coefficients of exactly one
% scaling function. All, except the `vm` first and last columns, contains
% translated versions of the same scaling function. The returned matrix have
% size [N, M]. Thus this function is ideal if you need to plot a wavelet
% approximation to a function, given some wavelet coefficents
%
% INPUT
% vm - Number of vanishing moments for a Daubechies wavelet
% N - Number of row in output the matrix
% M - Number of columns in the output matrix
% It is assumed that M = 2^r and N = 2^(r+q) for positive integers r and q
%
% RETURNS
% A N × M matrix
function A = cww_get_scaling_matrix_bd(log2N, log2M, wname)

    N = 2^log2N;
    M = 2^log2M;

    nres = log2N - log2M;
    vm = cww_extract_vm_from_wname(wname);

    % Generate scaling function
    j0 = cww_compute_j0(vm);
    r = j0 + nres;
    N_inner = 2^r;
    a = 2^nres;

    phi = zeros([N_inner,1]);
    phi(vm+1) = 2^(nres/2);
    phi = wl_idwt_impl(phi, wname, 'm',  nres, 'bd_mode', 'per' );
    phi = phi';
    ub = a*(2*vm);
    phi = phi(a+1:ub);

    x = zeros([N,1]);
    x(1:vm*a) = phi(a*(vm-1)+1:a*(2*vm-1));
    x(N-a*(vm-1)+1:N) = phi(1:a*(vm-1));

    A = zeros([N,M]);

    % Generate the inner matrix first
    for i = vm+1:M-vm
        A(:,i) = circshift(x, a*(i-1));
    end

    % Generate left edge boundary functions
    for k = 1:vm
        phi = zeros([N_inner,1]);
        phi(k) = 2^(nres/2);
        phi = wl_idwt_impl(phi, wname, 'm', nres, 'bd_mode', 'bd', 'prefilter_mode', 'bd_pre');
        phi = phi;
        ub = a*(vm+k-1);
        phi = phi(1:ub);
        A(1:ub, k) = phi;
    end

    % Generate right edge boundary functions
    for k = 1:vm
        phi = zeros([N_inner,1]);
        phi(2^(j0)-k+1) = 2^(nres/2);
        phi = wl_idwt_impl(phi, wname, 'm', nres, 'bd_mode', 'bd', 'prefilter_mode', 'bd_pre');
        phi = phi;
        phi = phi(N_inner-a*(vm+k-1)+1:N_inner);
        A(N-a*(vm+k-1)+1:N, M-k+1) = phi;
    end

    A = 2^(log2M/2)*A;
end

% Creates a matrix where each column contain the coefficients of exactly one
% scaling function. All columns contain translated versions of the same scaling 
% function. The returned matrix have size [N, M]. Thus this function is ideal 
% if you need to plot a wavelet approximation to a function, given some wavelet 
% coefficents
%
% INPUT
% vm - Number of vanishing moments for a Daubechies wavelet
% N - Number of row in output the matrix
% M - Number of columns in the output matrix
% It is assumed that M = 2^r and N = 2^(r+q) for positive integers r and q
%
% RETURNS
% A N × M matrix
function A = cww_get_scaling_matrix_per(log2N, log2M, wname)

    dwtmode('per', 'nodisp');

    N = 2^log2N;
    M = 2^log2M;

    nres = log2N - log2M;
    vm = cww_extract_vm_from_wname(wname);

    j0 = cww_compute_j0(vm);
    r = j0 + nres;
    N_inner = 2^r;
    a = 2^nres;

    %phi = zeros([N_inner,1]);
    %phi(vm+1) = 2^(nres/2);
    %phi = idwt_impl(phi, wname, nres, 'per');
    %phi = phi';
    %ub = a*(2*vm);
    %phi = phi(a+1:ub);

    S = csl_get_wavedec_s(r,nres);
    psi = zeros([1,N_inner]);
    psi(vm) = 2^(nres/2);
    psi = waverec(psi, S, wname);
    ub = 2^(nres)*(2*vm - 1);
    psi = psi(1:ub);
    psi = 2^(log2M/2)*psi;

    x = zeros([N,1]);
    x(1:2^nres * vm) = psi(2^nres*(vm-1)+1:2^nres*(2*vm-1));
    x(N-2^nres*(vm-1)+1:N) = psi(1:2^nres*(vm-1));
    A = zeros([N,M]);

    for i = 1:M
        A(:,i) = circshift(x, 2^nres*(i-1));
    end
    %A = 2^(log2M/2)*A;
end










