%% used to compute size of psatz
na_g = 3;
nb_g = 2;
na_c = 4;
nb_c = 3;
N = na_g + nb_g;
T = 10;
% Full method
d = 2;
simga_f = nchoosek(2*(N+T)-1+d,d)
zeta_f = nchoosek(2*(N+T)-1+d-1,d-1)
psi_f = zeta_f
mu_f = nchoosek(2*(N+T)-1+2*d-2,2*d-2)

% Alternative method
d = 1;
sigma_a = nchoosek(N+d,d)
zeta_a = sigma_a
psi_a = zeta_a
mu_a = nchoosek(N+2*d-1,2*d-1)

