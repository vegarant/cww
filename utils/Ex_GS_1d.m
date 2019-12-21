clear all
R = 5;
q = 1;
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

f = @(x) 50*(x-0.5).*(x-0.05).*(x-0.95) + x;%cos(2*pi*x)  + 0.2 * cos(2*pi*8 *x)+0.5*x; 
%f = @(x) (x+1).*(x-0.5).*(x-0.25).*(x+3).*(x-0.6); %+ cos(2*pi*x).*(x <= 0.5); 

samples = cww_walsh_sampling_1d(f,N);

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);

G = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 

A = cww_get_scaling_matrix(R+q+3, R, wname, bd_mode);
t = linspace(0,1,8*N);

sc = lsqr(G,samples, [], 600);

x = A*sc; 

figure()
eps = 1e-14;
t1 = linspace(0,1-eps,8*N)';
plot(t1,f(t1));
hold('on');
%t = linspace(0,1-eps,N)';
plot(t,x);
xlabel('t');
title(sprintf('vm: %d, N: %d, M: %d', vm, N, M));

% Truncated Walsh series
s = 0;
for n = 1:N
    s = s + samples(n)*wal(n-1, t1);
end

hold('on');
plot(t1,s);

%
%%legend({'f(t)', 'GS rec.'});
%%legend({'f(t)', 'GS rec.', 'f_trunk'});
%legend({'f(t)', 'f_trunk'});
%
%fprintf('||f-f_{N,M}||_{inf}: %g\n', norm(f(t)-x, inf));
%fprintf('||f-f_trunk||_{inf}: %g\n', norm(f(t)-s', inf));











