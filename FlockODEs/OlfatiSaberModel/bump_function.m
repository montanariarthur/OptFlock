function rhoh = bump_function(z,param)
% Equation (10) of Reference [2].
% h must lie in (0,1)

h = param.h;

if z >= 0 && z < h
    rhoh = 1;
elseif z >= h && z <= 1
    rhoh = 0.5*(1 + cos( pi*(z - h)/(1 - h) ));
else
    rhoh = 0;
end

end