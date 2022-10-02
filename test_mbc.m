%% test model-based control
clear
clc
yalmip('clear')
rng(1)

% generate system
% system in z, should be strictly proper
z = tf('z',0.1);     % sampling time 0.1
% Gz = z/(1.1-z)/(1.1+z);
Gz = z/(1.1-z)/(1.1+z)/(z+0.5);

% transform to system in lambda
[Gl,al,bl] = sys_trans(Gz);     % al, bl (low to high order)
na_g = length(al);              % size of model a
nb_g = length(bl);              % size of model b
na_c = 4;                       % size of controller a
nb_c = 3;                       % size of controller b

% define closed-looo acl = (1+A)(1+At)+BBt
ac = sdpvar(na_c,1);
bc = sdpvar(nb_c,1);
da = na_g+na_c;
db = nb_g+nb_c;
AA = compute_coeff([1;al],[1;ac]);
BB = compute_coeff([0;bl],[0;bc]);
AA(1) = [];
BB(1) = [];
if length(AA) > length(BB)
    BB = [BB; zeros(da-db,1)];
end
acl = AA + BB;

% solve MBC
opts = sdpsettings('solver','mosek','verbose', 0);
sol = optimize([], norm(acl, 1), opts)

% extract solution
gamma = value(norm(acl, 1))
ac = value(ac);
bc = value(bc);

[Cz, ACL_z, Cl, ACL_l] = recover_sol(Gl, ac, bc)
[zero, pole] = zpkdata(ACL_z);
zero = zero{1};
pole = pole{1}

