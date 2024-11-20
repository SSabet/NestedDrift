function noarb = FOC_bdotzero_resid(b,a,z,d,Va,forward,par)
    
    % Determine consumption from the passed information
    c = driftLiquid(0,d,b,a,z,par);

    % And evaluate the no-arbitrage condition with these policies
    noarb = MU(c,par).*(1 + par.chi0 * (2*forward - 1) + par.chi1 * d /a) - Va;

end