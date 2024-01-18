function eq = two_asset_kinked_cost(d,a, chi0, chi1)
    eq = chi0.*abs(d) + chi1.*d.^2/2.*(max(a,1e-6)).^(-1);
end