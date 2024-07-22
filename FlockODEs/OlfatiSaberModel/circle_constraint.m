function [c, ceq] = circle_constraint(x, c, r)
    % Circle constraint
    c = (x(1) - c(1))^2 + (x(2) - c(2))^2 - r^2;
    ceq = [];
end