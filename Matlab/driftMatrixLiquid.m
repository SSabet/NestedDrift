function BB = driftMatrixLiquid(s, dbb, dbf, par)
    
    % Unpack parameter values
    I = par.I; J = par.J; Nz = par.Nz;
    
    % Prepare elements for diag inputs
    X = -min(s,0)./dbb;
    Y =  min(s,0)./dbb - max(s,0)./dbf;
    Z =  max(s,0)./dbf;
    
    % Build the sparse drift matrix information
    rows = []; cols = []; vals = [];
    for k = 1:Nz
       for j = 1:J
           base = (I*J*(k-1))+(I*(j-1));
           rows = [rows, base+(1:I), base+(1:(I-1)), base+(2:I)]; % Order is central diag, updiag, lowdiag
           cols = [cols, base+(1:I), base+(2:I),     base+(1:(I-1))];
           vals = [vals; squeeze(Y(:,j,k)); squeeze(Z(1:(I-1),j,k)); squeeze(X(2:I,j,k))];
       end
    end
    
    % ...and make a sparse matrix out of it
    BB = sparse(rows,cols,vals,I*J*Nz,I*J*Nz);
    
return