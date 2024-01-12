function [a,daf,dab] = exponentialGrid(amin,amax,anum,coeff,power)

    % Build nonlinear grid
    x    = linspace(0,1,anum);
    xx   = x + coeff*x.^power;
    a    = amin + (amax-amin)/(max(xx) - min(xx))*xx;
    
    % Build difference grids
    daf = ones(anum,1); dab = ones(anum,1);
    daf(1:end-1) = diff(a);
    daf(end)     = daf(end-1); 
    
    dab(2:end)   = diff(a);
    dab(1)       = dab(2);


end