function noarb = fone(b,a,z,d,Va,forward,chi0,chi1,ga,xi,w,rb)

   noarb = (rb*b + (1-xi)*w*z - d - two_asset_kinked_cost(d,a))^(-ga) * (1 + chi0 * (2*forward - 1) + chi1 * d ^(chi1-1)/a) - Va;
    
end