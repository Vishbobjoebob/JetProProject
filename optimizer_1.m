function [x_opt, fval] = optimizer_1(options, x0, LB, UB, prm)

    f_obj_1 = @(X)my_obj_1(X, prm);
    f_con_1 = @(X)my_con_1(X, prm);

    [x_opt, fval] = fmincon(f_obj_1, x0, [], [], [], [], LB, UB, f_con_1, options);

end

function F_obj = my_obj_1(X, prm)

    
    [TSFC, ~, ~, ~, ~, ~, ~, ~] = engine_outputs(prm(1:3), X(1), X(2), X(3), X(4), X(5), X(6), X(7), prm(4));

    F_obj = TSFC;
end

function [Cneq, Ceq] = my_con_1(X, prm)

    [~, ST, f_max_main, f_max_ib, f_max_ab, ~, ~, ~] = engine_outputs(prm(1:3), X(1), X(2), X(3), X(4), X(5), X(6), X(7), prm(4));
    
    % Compressor Pressure Ratio
    Cneq(1) = X(1) - (55/X(2));

    % Maximum Fuel
    Cneq(2) = -f_max_main + X(5);
    Cneq(3) = -f_max_ib + X(6);
    Cneq(4) = -f_max_ab + X(7);
    
    % Specific Thrust
    Ceq(1) = ST - prm(5);
    
end


