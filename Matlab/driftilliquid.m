function sa = driftilliquid(d,a,z,par)

    % Calculate drift
    sa = par.ra.*a + par.xi.*par.w.*z + d;

end