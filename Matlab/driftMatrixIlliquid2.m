function KK = driftMatrixIlliquid2(sa,aaa,par)
    
    % Unpack parameter values
    I = par.I; J = par.J; Nz = par.Nz;
    
    % preallocating
    X = zeros(I,J,Nz);
    Z = zeros(I,J,Nz);

    % Prepare elements for diag inputs
    X(:,2:J,:) = -min(sa(:,2:J,:),0)./(aaa(:,2:J,:)-aaa(:,1:J-1,:));
    Z(:,1:J-1,:) =  max(sa(:,1:J-1,:),0)./(aaa(:,2:J,:)-aaa(:,1:J-1,:));
    Y = -(X+Z);

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