% Implementation of 
% "Superstabilizing Control of Discrete-Time ARX Models under Error in Variables"
% Full method
% cannot deal with large size due to complexity

clear
yalmip('clear')
rng(1)
tic

%% parameter
eps = [0.0;0.0];            % noise bound    dy,du
T = 4;                      % # of samples
d = 2;                      % degree of psatz
opts = sdpsettings('solver','mosek','verbose', 0);

%% generate system
% system in z, should be strictly proper
z = tf('z',0.1);     % sampling time 0.1
Gz = z/(1.1-z)/(1.1+z);

%% transform to system in lambda
[Gl,al,bl] = sys_trans(Gz);    % al, bl (low to high order, follows paper)
na_g = length(al);
nb_g = length(bl);
na_c = 1;                       % size of controller a
nb_c = 1;                       % size of controller b
n = [na_g;nb_g;na_c;nb_c];

%% generate trajectory
y = zeros(T+na_g,1);
y(1:na_g) = rand(na_g,1); 
u = rand(T+nb_g-1,1);        

for i = 1:T
    y(i+na_g) = -al'*y(i+na_g-1:-1:i) + bl'*u(i+nb_g-1:-1:i);
end

y_noise = y + eps(1)*(2*rand(T+na_g,1));
u_noise = u + eps(2)*(2*rand(T+nb_g-1,1));

sim = struct('y_noise',y_noise,'u_noise',u_noise,'epsilon',eps);

%% solve SS
out = Full_SS(sim, n, d, T);
sol = optimize(out.cons, out.obj, opts)

%% extract solution
gamma = value(out.obj)
ac = value(out.ac);
bc = value(out.bc);

[Cz, ACL_z, Cl, ACL_l] = recover_sol(Gl, ac, bc)
t = toc
