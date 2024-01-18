function eq = two_asset_kinked_cost_s(d,a, chi0, chi1)
    eq = zeros(size(a));
    eq(a ~= 0) = chi0.*abs(d(a ~= 0)) + chi1.*d(a ~= 0).^2.*(1./a(a~=0))/2;
end