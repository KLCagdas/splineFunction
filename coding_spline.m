val_FN = 3.2;
zt = 1:1:5;
syms a1 a2 a3 a4 a5 b1 b2 b3 b4 b5 c1 c2 c3 c4 c5;
a = [a1 a2 a3 a4 a5];
b = [b1 b2 b3 b4 b5];
c = [c1 c2 c3 c4 c5];

ab = [a; b; zeros(lenzt, 1)];
for i=1:lenzt-1
  splineDer = ab(i) - ab(i+1);

lenzt = length(zt);
data = [zt', zeros(lenzt, 1)];
Fzt1 = [a .* zt.^2; b.* zt; c]';
Fzt2 = zeros(5, 3);

for i=2:lenzt
  Fzt2(i-1) += a(i) * zt(i-1)^2 + b(i) * zt(i-1) + c(i);
end

function find_zt 
# loop for storing the values of FN in respect to zt
for i=1:lenzt
  dx = x(n) - xt(n);    dz = zi(n) - zt(n);
  rho = ((x(n) - xt(n))^2 + (z(n) - zt(n))^2);
  cEq = 24*E*s^6;
  deriv_LenJones = cEq * (rho^-4 - 2 * s^6 * rho^-7);
  FN = ((z(i) - zt(i)) * deriv_LenJones);
  data(:, i) = FN;
end
# now data is a matrix of zt, FN 

  Fzt1 = [Fzt1, FN];
end