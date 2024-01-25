function KK = driftMatrixIlliquid(sa,aaa,par)
    
    % Unpack parameter values
    I = par.I; J = par.J; Nz = par.Nz;
    diag = 1:(I*J*Nz); diagArray = reshape(diag,I,J,Nz);
    
    % preallocating
    X = zeros(I,J,Nz);
    Z = zeros(I,J,Nz);

    % Prepare elements for diag inputs
    X(:,2:J,:)   = -min(sa(:,2:J,:),0)./(aaa(:,2:J,:)-aaa(:,1:J-1,:));
    Z(:,1:J-1,:) =  max(sa(:,1:J-1,:),0)./(aaa(:,2:J,:)-aaa(:,1:J-1,:));
    Y            = -(X+Z);
    
    % Build components (order is diag, upper, lower)
    row = [reshape(diagArray(:,(2:end),:),1,[])   diag            reshape(diagArray(:,1:(end-1),:),1,[])  ]; 
    col = [reshape(diagArray(:,(2:end),:),1,[])-I diag            reshape(diagArray(:,1:(end-1),:),1,[])+I];
    val = [reshape(X(:,2:end,:),1,[])             reshape(Y,1,[]) reshape(Z(:,1:(end-1),:),1,[])          ];
    
    % Make sparse transition matrix
    KK = sparse(row,col,val,I*J*Nz,I*J*Nz);
    
return
