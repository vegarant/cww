% This script compares generalised sampling and truncated walsh reconstruction in 
% two dimensions, using a finite number of samples. By adjusting the number of 
% vanishing momenents one can see how the reconstruction improves.

clear all

load('cww_defaults.mat') % load font size, line width, etc.


R = 4;
q = 1;
q2 = 6; % Degree of high resolution image
vms = {2,4,6};
is_per = 0;
dest = 'plots';
do_save = 1;

M = 2^R;
N = 2^(R+q);
log2N = R+q;
log2M = R;

if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

N2 = 2^(R+q+q2)
t2 = linspace(0,1,N2);

[X2,Y2] = meshgrid(t2,t2);

%f = @(x,y) 5*(x-1).*(x-0.5).*(x-0.25).*cos(4*pi*y) + 5;
f = @(x,y) cos(1.5*pi*x).*sin(3*pi*y) ;%+ (x>= 0.5);
F2 = f(X2,Y2);

samples = cww_sample_walsh_2d(f, N);

F_approx = zeros(N2,N2);
F_approx(1:N, 1:N) = samples;
F_approx = cww_fastwht_2d(F_approx)*N2*N2;

max_diff = max(abs(F_approx(:)-F2(:)))
n = length(vms);
W_approxs = cell(n,1);

for i = 1:length(vms)

    vm = vms{i};
    wname = sprintf('db%d', vm);
    j0 = cww_compute_j0(vm);

    phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
    G = @(x,mode) cww_handle_2d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 

    sc = lsqr(G, samples(:), 1e-10, 600);
    sc = reshape(sc, M, M);

    W = cww_map_wcoeff_to_func_vals_2d(sc, R+q+q2, wname, bd_mode);
    W_approxs{i} = W;
    mdiff = max(abs(W(:) - F2(:)));
    if mdiff > max_diff
        max_diff = mdiff;
    end
end

n_F = norm(F2, 'fro');
n_F_app = norm(F2 - F_approx, 'fro');

mi = min(F2(:));
ma = max(F2(:));
F_plot = (F2 - mi)/(ma-mi);
F_approx_plot = (F_approx - mi)/(ma-mi);


for i = 1:n
    vm = vms{i};
    wname = sprintf('db%d', vm);
    W = W_approxs{i};
    W_plot = (W - mi)/(ma-mi);
    
    % Compute error 
    n_W_app = norm(F2 - W, 'fro');
    fprintf('TW: rel err: %.4f, Wavelet (%s), rel err: %.4f\n', n_F_app/n_F, wname, n_W_app/n_F);
    W_error = 1-(abs(F2-W)/max_diff);
    
    fname_W = sprintf('comp_2d_GS_N_%d_M_%d_%s.%s', N, M, wname, cww_dflt.image_format);
    fname_W_error = sprintf('comp_2d_GS_N_%d_M_%d_%s_error.%s', N, M, wname, cww_dflt.image_format);

    if do_save
    
        imwrite(im2uint8(W_plot), fullfile(dest, fname_W));
        imwrite(im2uint8(W_error), fullfile(dest, fname_W_error));
    
    end
    
end

fname_func = sprintf('comp_2d_smooth_org_func.%s', cww_dflt.image_format);
fname_F_approx = sprintf('comp_2d_TW_N_%d.%s', N, cww_dflt.image_format);
fname_F_approx_error = sprintf('comp_2d_TW_N_%d_error.%s', N, cww_dflt.image_format);

F_approx_error = 1-(abs(F2-F_approx)/max_diff);

if do_save

        imwrite(im2uint8(F_plot), fullfile(dest, fname_func));
        imwrite(im2uint8(F_approx_plot), fullfile(dest, fname_F_approx));
        imwrite(im2uint8(F_approx_error), fullfile(dest, fname_F_approx_error));

end



