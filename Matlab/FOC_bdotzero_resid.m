function noarb = FOC_bdotzero_resid(b,a,zk,d,Va,forward,rb, MU,par)

   noarb = MU(rb*b + (1-par.xi)*par.w*zk - d - two_asset_kinked_cost(d,a, par.chi0, par.chi1)).*...
        (1 + par.chi0 * (2*forward - 1) + par.chi1 * d /a) - Va;
end