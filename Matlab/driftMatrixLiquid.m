function BB = driftMatrixLiquid(sb, bbb, par)
    
    % Unpack parameter values
    I = par.I; J = par.J; Nz = par.Nz;
    diag = 1:(I*J*Nz);
    
    % Preallocating memory
    X = zeros(I,J,Nz);
    Z = zeros(I,J,Nz);

    % Prepare elements for diag inputs
    X(2:I,:,:)   = -min(sb(2:I,:,:),0)  ./(bbb(2:I,:,:)- bbb(1:I-1,:,:));
    Z(1:I-1,:,:) =  max(sb(1:I-1,:,:),0)./(bbb(2:I,:,:)- bbb(1:I-1,:,:));
    Y            = -(X+Z);
    
    % Reshape inputs to long
    longX = reshape(X,1,[]);
    longZ = reshape(Z,1,[]);
    
    % Build elements
    rows = [diag(2:end)     diag            diag(1:(end-1))];
    cols = [diag(2:end)-1   diag            diag(1:(end-1))+1];
    vals = [longX(2:end)    reshape(Y,1,[]) longZ(1:(end-1))];
    
    % Build sparse transition matrix
    BB = sparse(rows,cols,vals,I*J*Nz,I*J*Nz);
    
return
