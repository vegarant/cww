clear all

load('cww_defaults.mat') % load font size, line width, etc.


R = 4;
q = 1;
q2 = 3; % Degree of high resolution image
vm = 2;
is_per = 0;
dest = 'plots';
do_save = 1;

M = 2^R;
N = 2^(R+q);
wname = sprintf('db%d', vm);
j0 = cww_compute_j0(vm);
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
phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
G = @(x,mode) cww_handle_2d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 

sc = lsqr(G, samples(:), [], 600);
sc = reshape(sc, M, M);

W = cww_map_wcoeff_to_func_vals_2d(sc, R+q+q2, wname, bd_mode);

F_approx = zeros(N2,N2);
F_approx(1:N, 1:N) = samples;
F_approx = cww_fastwht_2d(F_approx)*N2*N2;



%figure();
%subplot(131); imagesc(F2); colormap('gray'); colorbar(); title('Original');
%subplot(132); imagesc(W); colormap('gray'); colorbar(); title('Wavelet approximation');
%subplot(133); imagesc(F_approx); colormap('gray'); colorbar(); title('Walsh approximation');

mi = min(F2(:));
ma = max(F2(:));

F_plot = (F2 - mi)/(ma-mi);
F_approx_plot = (F_approx - mi)/(ma-mi);
W_plot = (W - mi)/(ma-mi);

F_approx_plot(F_approx_plot > 1) = 1;
F_approx_plot(F_approx_plot < 0) = 0;

W_plot(W_plot > 1) = 1;
W_plot(W_plot < 0) = 0;

fname_func = sprintf('comp_2d_smooth_org_func.%s', cww_dflt.image_format);
fname_F_approx = sprintf('comp_2d_TW_N_%d.%s',N, cww_dflt.image_format);
fname_W = sprintf('comp_2d_GS_N_%d_M_%d_%s.%s', N, M, wname, cww_dflt.image_format);

if do_save

    imwrite(im2uint8(F_plot), fullfile(dest, fname_func));
    imwrite(im2uint8(F_approx_plot), fullfile(dest, fname_F_approx));
    imwrite(im2uint8(W_plot), fullfile(dest, fname_W));

end







