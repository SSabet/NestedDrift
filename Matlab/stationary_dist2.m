function g = stationary_dist2(A, a, b, par)
    % Unpack parameter values
    I = par.I; J = par.J; Nz = par.Nz;
    
    M          = I*J*Nz;
    AT         = A';
    
    % building a matrix of da*db*dz for aggregation
    db = ones(I,1);
    da = ones(J,1);
    
    db(2:I-1) = (b(3:I) - b(1:I-2))/2;
    db(1) = b(2)-b(1);
    db(I) = b(I) - b(I-1);

    da(2:J-1) = (a(3:J) - a(1:J-2))/2;
    da(1) = a(2) - a(1);
    da(J) = a(J) - a(J-1);
    
    dadb = db*da';
    dadbdz = repmat(dadb, 1,1,Nz);

    % Fix one value so matrix isn't singular:
    vec        = zeros(M,1);
    iFix       = 10;
    vec(iFix)  = .01;
    AT(iFix,:) = [zeros(1,iFix-1),1,zeros(1,M-iFix)];
    
    % Solve system:
    g_stacked  = AT\vec;
    g          = reshape(g_stacked,I,J,Nz); % reshape the density

    % Normalise total distibution to 1
    g_sum = sum(sum(sum(g.*dadbdz)));
    % g_sum      = g_stacked'*ones(M,1)*da*db;
    % g_stacked  = g_stacked./g_sum;
    g = g/g_sum;
        
end