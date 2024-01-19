function noarb = FOC_bdotzero_resid(b,a,zk,d,Va,forward,rb,par)

   noarb = MU(rb*b + (1-par.xi)*par.w*zk - d - adjustmentCost(d,a, par.chi0, par.chi1),par).*...
        (1 + par.chi0 * (2*forward - 1) + par.chi1 * d /a) - Va;
end