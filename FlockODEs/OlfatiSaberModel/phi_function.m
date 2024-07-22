function phi = phi_function(z,param)
% Equation (15) of Reference [2].

% 0 < a <= b
a = param.a;
b = param.b;
c = abs(a - b)/sqrt(4*a*b);

arg = z + c;
sigma1 = arg / sqrt(1 + arg^2);
phi = 0.5*( (a + b)*sigma1  + (a - b) );

end