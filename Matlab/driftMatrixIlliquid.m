function KK = driftMatrixIlliquid(m,dkb,dkf,par)
    
    % Unpack parameter values
    I = par.I; J = par.J; Nz = par.Nz;
    
    % Prepare elements for diag inputs
    X = -min(m,0)./dkb;
    Y =  min(m,0)./dkb - max(m,0)./dkf;
    Z =  max(m,0)./dkf;
    
    % Build the sparse drift matrix information
    row = []; col = []; val = [];
    for k = 1:Nz
            
        diag = I*J*(k-1)+(1:I);
        row  = [row, diag, diag];
        col  = [col, diag, diag+I];
        val  = [val; Y(:,1,k); Z(:,1,k)];
        
        for j=2:(J-1)
            diag = I*J*(k-1) + I*(j-1)+(1:I);
            row  = [row, diag, diag, diag];
            col  = [col, diag, diag+I, diag-I];
            val  = [val; Y(:,j,k); Z(:,j,k); X(:,j,k)];
        end
        
        diag = I*J*(k-1) + I*(J-1)+(1:I);
        row  = [row, diag, diag];
        col  = [col, diag, diag-I];
        val  = [val; Y(:,J,k); X(:,J,k)];
    end
    
    % ...and make a sparse matrix out of it
    KK = sparse(row,col,val,I*J*Nz,I*J*Nz);
    
return