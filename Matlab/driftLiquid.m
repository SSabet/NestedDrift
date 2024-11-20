function sb = driftLiquid(c,d,b,a,z,par)

    % Calculate state-dependent rate of return
    Rb = par.rb_pos.*(b>0) + par.rb_neg.*(b<0);

    % Calculate drift
    sb = (1-par.xi).*par.w.*z + Rb.*b - d - adjustmentCost(d,a,par.chi0,par.chi1) - c;

end