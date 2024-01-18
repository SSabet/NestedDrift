function u = U(c,par)

    gamma = par.gamma;
    
    % Utility function, marginal utility and its inverse
    if gamma==1
        u = log(c);
    else
        u = c.^(1-gamma)/(1-gamma);
    end

end