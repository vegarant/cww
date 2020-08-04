% This script is not working.

clear all
R = 4;
q = 1;
q2 = 5;
vm = 4;
is_per = 0;

M = 2^R;
N = 2^(R+q);
N2 = 2^(R+q+q2);
wname = sprintf('db%d', vm);
j0 = cww_compute_j0(vm);
log2N = R+q;
log2M = R;

if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

f = @(x) 50*(x-0.6).*(x-0.5).*(x-0.5).*(x-0.05).*(x-0.95) + x + 0.5*(x > 0.5-(1/(2*N)));%cos(2*pi*x)  + 0.2 * cos(2*pi*8 *x)+0.5*x; 
%f = @(x) (x+1).*(x-0.5).*(x-0.25).*(x+3).*(x-0.6); %+ cos(2*pi*x).*(x <= 0.5); 

r = 5;
samples_cont = cww_walsh_sampling_1d(f,N, r);
r = 0;
samples_disc = cww_walsh_sampling_1d(f,N, r);

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);

G = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 

A = cww_get_scaling_matrix(R+q+q2, R, wname, bd_mode);
t = linspace(0,1,N2);

sc_cont = lsqr(G,samples_cont, [], 600);
sc_disc = lsqr(G,samples_disc, [], 600);

x_cont = A*sc_cont; 
x_disc = A*sc_disc; 

figure()
eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
plot(t1,f(t1));
hold('on');
%t = linspace(0,1-eps,N)';
plot(t,x_cont);
xlabel('t');
title(sprintf('vm: %d, N: %d, M: %d', vm, N, M));

%hold('on');
%%t = linspace(0,1-eps,N)';
%plot(t,x_disc);
%xlabel('t');

raw_samples = zeros([N2,1]);
raw_samples(1:N) = samples_cont;
walsh_approx = fastwht(raw_samples)*N2;
plot(t1, walsh_approx)

%hold('on')
%raw_samples(1:N) = samples_disc;
%walsh_approx = fastwht(raw_samples)*N2;
%plot(t1, walsh_approx)

legend({'f(t)', 'GS\_cont', 'wal\_cont' })
%legend({'f(t)', 'GS\_cont', 'GS\_disc', 'wal\_cont', 'wal\_disc'})
