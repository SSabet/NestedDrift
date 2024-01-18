function g = stationary_dist(A, a, b, par)

    % Unpack parameter values
    I = par.I; J = par.J; Nz = par.Nz;
    
    M          = I*J*Nz;
    AT         = A';
    
%     % building a matrix of da*db*dz for aggregation
%     db = ones(I,1);
%     da = ones(J,1);
%     
%     db(2:I-1) = (b(3:I) - b(1:I-2))/2;
%     db(1) = b(2)-b(1);
%     db(I) = b(I) - b(I-1);
% 
%     da(2:J-1) = (a(3:J) - a(1:J-2))/2;
%     da(1) = a(2) - a(1);
%     da(J) = a(J) - a(J-1);
%     
%     dadb = db * da';
%     dadbdz = repmat(dadb, 1,1,Nz);

    % Fix one value so matrix isn't singular:
    vec        = zeros(M,1);
    iFix       = 10;
    vec(iFix)  = .01;
    AT(iFix,:) = [zeros(1,iFix-1),1,zeros(1,M-iFix)];
    
%     % Solve system:
%     g_stacked  = AT\vec;
%     g          = reshape(g_stacked,I,J,Nz); % reshape the density
% 
%     % NOTE FROM PATRICK: is g_stacked the density or the historgram? 
%     
%     % Normalise total distibution to 1
%     g_sum = sum(sum(sum(g.*dadbdz)));
%     % g_sum      = g_stacked'*ones(M,1)*da*db;
%     % g_stacked  = g_stacked./g_sum;
%     g = g/g_sum;
%     
%     assert(sum(sum(sum(g)))==1, "Normalisation error: this isn't actually a density")

    % Alternative approach (Patrick)
        
    % Solve system and rebase to get the **histogram** over the state space:
    gtilde_stacked = AT\vec; gtilde_stacked = gtilde_stacked./sum(gtilde_stacked);
    
    % Prepare objects to convert the histogram to density in case grids are
    % nonlinear. This uses a 2-dimensional equivalent of the method outlined in
    % Section 7 of the Appendix to Achdou et al HACT paper
    db_tilde = ones(I,1); 
    db_tilde(2:I-1) = (b(3:I) - b(1:I-2))/2; 
    db_tilde(1) = b(2)-b(1); 
    db_tilde(I) = b(I)-b(I-1);
    
    da_tilde = ones(J,1); 
    da_tilde(2:J-1) = (a(3:J) - a(1:J-2))/2; 
    da_tilde(1) = a(2)-a(1); 
    da_tilde(J) = a(J)-a(J-1);
    
    [AAA,BBB] = meshgrid(da_tilde,db_tilde); 
    dba_tilde = reshape(AAA .* BBB,I*J,1);
    grid_diag = spdiags(repmat(dba_tilde,Nz,1),0,M,M);
    
    % Convert the histogram to a density
    g_stacked = grid_diag \ gtilde_stacked;
    
    % Reshape to return
    g = reshape(g_stacked,I,J,Nz);
    
    assert(sum(sum(sum(g)))==1, "Normalisation error: this isn't actually a density")
        
end
