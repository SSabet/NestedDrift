function noarb = fone(b,a,zk,d,Va,forward,rb,par)
   cellfun(@(x) assignin('caller', x, par.(x)), fieldnames(par));
   noarb = (rb*b + (1-xi)*w*zk - d - two_asset_kinked_cost(d,a, chi0, chi1)).^(-par.gamma).*...
       (1 + chi0 * (2*forward - 1) + chi1 * d ^(chi1-1)/a) - Va;
    
end