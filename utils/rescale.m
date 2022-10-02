function sys2 = rescale(sys)
% rescale such that the coefficient of highest order term = 1
[b,a] = tfdata(sys);
b = b{1};
a = a{1};
b = b/a(1);
a = a/a(1);
sys2 = tf(b,a,0.1);
end