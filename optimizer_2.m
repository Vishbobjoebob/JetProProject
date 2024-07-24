function [y_opt, fval] = optimizer_2(options, x0, LB, UB, prm, x_opt)

    f_obj_2 = @(X)my_obj_2(X, prm, x_opt);
    f_con_2 = @(X)my_con_2(X, prm, x_opt);

    [y_opt, fval] = fmincon(f_obj_2, x0, [], [], [], [], LB, UB, f_con_2, options);

end

function F_obj = my_obj_2(X, prm, x_opt)

    
    [~, ST, ~, ~, ~, ~, ~, ~] = engine_outputs(prm(1:3), x_opt(1), x_opt(2), x_opt(3), X(1), X(2), X(3), X(4), prm(4));

    F_obj = -ST;
end

function [Cneq, Ceq] = my_con_2(X, prm, x_opt)

    [~, ~, f_max_main, f_max_ib, f_max_ab, ~, ~, ~] = engine_outputs(prm(1:3), x_opt(1), x_opt(2), x_opt(3), X(1), X(2), X(3), X(4), prm(4));
    
    % Maximum Fuel
    Cneq(2) = -f_max_main + X(2);
    Cneq(3) = -f_max_ib + X(3);
    Cneq(4) = -f_max_ab + X(4);

    Ceq(1) = 0;
    
    
end
