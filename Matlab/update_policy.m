function [c,d] = update_policy(Va, Vb, MU_inv, a, d_bar, par)
    cellfun(@(x) assignin('caller', x, par.(x)), fieldnames(par));

    c = MU_inv(Vb);
    d = two_asset_kinked_FOC(Va,Vb,a, chi0, chi1);

end