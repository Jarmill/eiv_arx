function out = Full_SS(sim, n, d, T)
%% Fual method, Algorithm 1
%  sim:     sampled trajectories
%    n:     size of coefficients
%    d:     degree of the psatz
%    T:     # of samples for design

% extract data and size
na_g = n(1);
nb_g = n(2);
na_c = n(3);
nb_c = n(4);
y_noise = sim.y_noise(1:T+na_g);        % noisy output
u_noise = sim.u_noise(1:T+nb_g-1);      % noisy input
eps = sim.epsilon;                      % noise bound epsilon
Cons = [];                              % save constraints

% define variables
ag = sdpvar(na_g,1);
bg = sdpvar(nb_g,1);
ac = sdpvar(na_c,1);
bc = sdpvar(nb_c,1);
dy = sdpvar(T+na_g,1);           % delta x in paper
du = sdpvar(T+nb_g-1,1);
vars = [ag;bg;dy;du];

% define constraints
g = [eps(1)-dy;eps(1)+dy;eps(2)-du;eps(2)+du];   
for i = 1:T                   
    ht = y_noise(i+na_g) + ag'*y_noise(i+na_g-1:-1:i) - bg'*u_noise(i+nb_g-1:-1:i);
    h(i,1) = ht - ag'*dy(i+na_g-1:-1:i) + bg'*du(i+nb_g-1:-1:i) - dy(i+na_g);
end
cons_data = struct('ineq', g, 'eq', h);    % constraints from data

% define m, acl
da = na_g+na_c;
db = nb_g+nb_c;
for i = 1:da        % eq.(15c)
    m(i,1) = polynomial([ag;bg],2*d);
end
AA = compute_coeff([1;ag],[1;ac]);
BB = compute_coeff([0;bg],[0;bc]);
AA(1) = [];
BB(1) = [];
if length(AA) > length(BB)
    BB = [BB; zeros(da-db,1)];
end
acl = AA + BB;

% define non-negative polynomial eq.(15de)
gamma = sdpvar;
p1 = m-acl;
p2 = m+acl;
p3 = gamma-sum(m);
P = [p1(:);p2(:);p3(:)];    % P contains M-A-BK >= 0, M+A+BK >= 0, 1-delta-sum(M,2) >= 0
l = length(P);

% get constraints
Gram = cell(l,1);                       % save Gram matrix
Coef = cell(l,1);                       % save coef of mu
for i = 1:l
    [p_psatz, cons_psatz, Gram{i,1}, Coef{i,1}] = Full_SS_psatz(P(i), cons_data, d, vars);
    Cons = [Cons; cons_psatz];          % constraints from psatz
end
out.poly = p_psatz;
out.cons = Cons;
out.Gram = Gram;
out.Coef = Coef;
out.obj = gamma;
out.ac = ac;
out.bc = bc;
end