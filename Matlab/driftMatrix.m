function AA = driftMatrix(sb, bbb, sa,aaa, par)
    
    % Unpack parameter values
    I = par.I; J = par.J; Nz = par.Nz;
    
    % preallocating
    X = zeros(I,J,Nz);
    Z = zeros(I,J,Nz);

    % Prepare elements for diag inputs of the liquid asset drift
    X(2:I,:,:) = -min(sb(2:I,:,:),0)./(bbb(2:I,:,:)- bbb(1:I-1,:,:));
    Z(1:I-1,:,:) =  max(sb(1:I-1,:,:),0)./(bbb(2:I,:,:)- bbb(1:I-1,:,:));
    Y = -(X+Z);

    % Build the row, column, value entries for sparse liquid drift matrix
    rows = []; cols = []; vals = [];
    for k = 1:Nz
       for j = 1:J
           base = (I*J*(k-1))+(I*(j-1));
           rows = [rows, base+(1:I), base+(1:(I-1)), base+(2:I)]; % Order is central diag, updiag, lowdiag
           cols = [cols, base+(1:I), base+(2:I),     base+(1:(I-1))];
           vals = [vals; squeeze(Y(:,j,k)); squeeze(Z(1:(I-1),j,k)); squeeze(X(2:I,j,k))];
       end
    end
    
    
    % Prepare elements for diag inputs
    X = zeros(I,J,Nz);
    Z = zeros(I,J,Nz);

    X(:,2:J,:) = -min(sa(:,2:J,:),0)./(aaa(:,2:J,:)-aaa(:,1:J-1,:));
    Z(:,1:J-1,:) =  max(sa(:,1:J-1,:),0)./(aaa(:,2:J,:)-aaa(:,1:J-1,:));
    Y = -(X+Z);
    
    % Add the row, column, value entries for sparse illiquid drift matrix
    for k = 1:Nz
            
        diag = I*J*(k-1)+(1:I);
        rows  = [rows, diag, diag];
        cols  = [cols, diag, diag+I];
        vals  = [vals; Y(:,1,k); Z(:,1,k)];
        
        for j=2:(J-1)
            diag = I*J*(k-1) + I*(j-1)+(1:I);
            rows  = [rows, diag, diag, diag];
            cols  = [cols, diag, diag+I, diag-I];
            vals  = [vals; Y(:,j,k); Z(:,j,k); X(:,j,k)];
        end
        
        diag = I*J*(k-1) + I*(J-1)+(1:I);
        rows  = [rows, diag, diag];
        cols  = [cols, diag, diag-I];
        vals  = [vals; Y(:,J,k); X(:,J,k)];
    end

    % ...and make a sparse matrix out of it
    AA = sparse(rows,cols,vals,I*J*Nz,I*J*Nz);
    
return