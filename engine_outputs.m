function [TSFC, ST, f_max_main, f_max_ib, f_max_ab, nth, np, no] = engine_outputs(ambient, Pr_c, Pr_f, beta, b, f, f_ib, f_ab, separated)

    engine = engineDesigner(ambient, Pr_c, Pr_f, beta, b, f, f_ib, f_ab);

    if separated
        TSFC = engine.TSFC_sn;
        ST = engine.ST_sn;
        nth = engine.nth_sn;
        np = engine.np_sn;
        no = engine.n0_sn;
    else
        TSFC = engine.TSFC_cn;
        ST = engine.ST_cn;
        nth = engine.nth_cn;
        np = engine.np_cn;
        no = engine.n0_cn;
    end
    
    f_max_main = engine.f_max_main;
    f_max_ib = engine.f_max_ib;
    f_max_ab = engine.f_max_ab;
end