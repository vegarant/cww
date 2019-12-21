clear all
R = 5;
q = 1;
q2 = 2; % Degree of high resolution image
vm = 4;
is_per = 0;

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
f = @(x,y) 10*cos(1.5*pi*x).*sin(3*pi*y) + 10*(x>= 0.5);
F2 = f(X2,Y2);

samples = cww_walsh_sampling_2d(f, N);
phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
G = @(x,mode) cww_handle_2d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 

sc = lsqr(G, samples(:), [], 600);
sc = reshape(sc, M, M);

W = cww_map_2d_wcoeff_to_func_vals(sc, R+q+q2, wname, bd_mode);

F_approx = zeros(N2,N2);
F_approx(1:N, 1:N) = samples;
F_approx = cww_fastwht_2d(F_approx)*N2*N2;

figure();
subplot(131); imagesc(F2); colormap('gray'); title('Original');
subplot(132); imagesc(W); colormap('gray'); title('Wavelet approximation');
subplot(133); imagesc(F_approx); colormap('gray'); title('Walsh approximation');













