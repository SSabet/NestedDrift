function BB = driftMatrixLiquid(s, dbb, dbf, par)
    
    cellfun(@(x) assignin('caller', x, par.(x)), fieldnames(par));

    % Preallocate some memory
    BBi      = cell(Nz,1);
    updiag   = zeros(I*J,Nz);
    lowdiag  = zeros(I*J,Nz);
    centdiag = zeros(I*J,Nz);
    
    % Prepare elements for diag inputs
    X = -min(s,0)./dbb;
    Y = min(s,0)./dbb - max(s,0)./dbf;
    Z = max(s,0)./dbf;

    for i = 1:Nz
        centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
    end

    lowdiag(1:I-1,:) = X(2:I,1,:);
    updiag(2:I,:)    = Z(1:I-1,1,:);
    for j = 2:J
        lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];
        updiag(1:j*I,:)  = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
    end

    for nz = 1:Nz
        BBi{nz}=spdiags(centdiag(:,nz),0,I*J,I*J)+spdiags([updiag(:,nz);0],1,I*J,I*J)+spdiags([lowdiag(:,nz);0],-1,I*J,I*J);
    end

    BB = [BBi{1}, sparse(I*J,I*J); sparse(I*J,I*J), BBi{2}];
    
return