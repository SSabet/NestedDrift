function KK = driftMatrixIlliquid(m,dkb,dkf,par)

    cellfun(@(x) assignin('caller', x, par.(x)), fieldnames(par));

    % Preallocate some memory
    AAi        = cell(Nz,1);
    AAupdiag   = zeros(I*J,Nz);
    AAlowdiag  = zeros(I*J,Nz);
    AAcentdiag = zeros(I*J,Nz);
    
    chi  = -min(m,0)./dkb;
    yy   = min(m,0)./dkb - max(m,0)./dkf;
    zeta = max(m,0)./dkf;

    for nz=1:Nz
        %This will be the upperdiagonal of the matrix AAi
        AAupdiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
        for j=1:J
            AAupdiag=[AAupdiag;zeta(:,j,nz)];
        end

        %This will be the center diagonal of the matrix AAi
        AAcentdiag= yy(:,1,nz);
        for j=2:J-1
            AAcentdiag=[AAcentdiag;yy(:,j,nz)];
        end
        AAcentdiag=[AAcentdiag;yy(:,J,nz)];

        %This will be the lower diagonal of the matrix AAi
        AAlowdiag=chi(:,2,nz);
        for j=3:J
            AAlowdiag=[AAlowdiag;chi(:,j,nz)];
        end

        % Add up the upper, center, and lower diagonal into a sparse matrix
        AAi{nz} = spdiags(AAcentdiag,0,I*J,I*J)+spdiags(AAlowdiag,-I,I*J,I*J)+spdiags(AAupdiag,I,I*J,I*J);

    end

    KK = [AAi{1}, sparse(I*J,I*J); sparse(I*J,I*J), AAi{2}];

return