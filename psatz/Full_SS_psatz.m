function [p_psatz, cons_psatz, Gram, Coef] = Full_SS_psatz(p, C, deg, vars)
%% Enforce Putinar Psatz using eq.(8)
%    p:   nonnegative polynomial
%    C:   eq and ineq constraints
%    d:   degree of p
% vars:   variables a,b

n_g = length(C.ineq);       % # of inequality
n_h = length(C.eq);         % # of equality
n_var = length(vars);       % # of variables
d_g = degree(C.ineq);       % degree of inequality
d_h = degree(C.eq);         % degree of equality
Gram = cell(1+n_g, 1);      % gram coefficients of sigma0, sigmai
Coef = cell(n_h, 1);        % coefficients of phi
p_psatz = 0;                % Putinar form of p

% generate sigma0 
[pow0, ~] = momentPowers(0, n_var, deg);   
vect_0 = recovermonoms(pow0, vars);     % % basis v(x)
l_0 = length(vect_0);
Gram{1} = sdpvar(l_0);                  % gram matrix of sigma0
sigma0 = vect_0'*Gram{1}*vect_0;        
p_psatz = p_psatz + sigma0;
cons_psatz = (Gram{1} >= 0):'Gram_0';   % psd constraint from sigma0

% generate sigmai
if mod(d_g,2) == 0                      % decide degree of poly
    [powi, ~] = momentPowers(0, n_var, floor((2*deg-d_g)/2));
else
    [powi, ~] = momentPowers(0, n_var, floor((2*deg-d_g-1)/2));
end
vect_i = recovermonoms(powi, vars);                         % basis v(x)
l_i = length(vect_i);
for i = 1:n_g      
    Gram{i+1} = sdpvar(l_i);                                % gram matrix of sigmai
    sigmai = vect_i'*Gram{i+1}*vect_i;
    p_psatz = p_psatz + sigmai*C.ineq(i); 
    cons_psatz = [cons_psatz; (Gram{i+1} >= 0):'Gram_i'];   % psd constraint from sigmai
end

% generate phi
[pow_phi, ~] = momentPowers(0, n_var, 2*deg-d_h);
vect_phi = recovermonoms(pow_phi, vars);            % basis v(x)
l_mu = length(vect_phi);    
for i = 1:n_h      
    Coef{i} = sdpvar(l_mu,1);                       % coefficients of phi
    phi = vect_phi'*Coef{i};
    p_psatz = p_psatz + phi*C.eq(i);   
end

cons_psatz = [cons_psatz; (coefficients(p - p_psatz, vars) == 0):'eq'];  % p = p_psatz
end