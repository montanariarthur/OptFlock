function qj = find_closest_point(qi, c, r)
% Equation (45) in Reference [2].

    % Define objective function
    objective_function = @(x) norm(x - qi);
    
    % Define constraints
    nonlcon = @(x) circle_constraint(x, c, r);

    % Initial guess
    x0 = [0; 0];

    % Solve using fmincon
    options = optimoptions('fmincon', 'Display', 'off');
    qj = fmincon(objective_function, x0, [], [], [], [], [], [], nonlcon, options);
end

