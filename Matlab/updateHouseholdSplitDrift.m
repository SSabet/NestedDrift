function sol = updateHouseholdSplitDrift(V,grids,par)

    % INCOMPLETE TIDYING. WOULD BE NICE TO MAKE THEM COMPARABLE THOUGH

    %Preallocation
    [VbF,VbB,VaF,VaB,c] = deal(zeros(par.I,par.J,par.Nz));
    
    updiag = zeros(I*J,Nz);
    lowdiag = zeros(I*J,Nz);
    centdiag = zeros(I*J,Nz);
    AAi = cell(Nz,1);
    BBi = cell(Nz,1);

    
    for n=1:maxit

        % DERIVATIVES W.R.T. b
        % forward difference
        VbF(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
        VbF(I,:,:)     = ((1-xi)*w*zzz(I,:,:) + Rb(I,:,:).*bmax).^(-ga); %state constraint boundary condition

        % backward difference
        VbB(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
        VbB(1,:,:) = ((1-xi)*w*zzz(1,:,:) + Rb(1,:,:).*bmin).^(-ga); %state constraint boundary condition
    
        %DERIVATIVES W.R.T. a
        % forward difference
        VaF(:,1:J-1,:) = (V(:,2:J,:)-V(:,1:J-1,:))/da;
        % backward difference
        VaB(:,2:J,:) = (V(:,2:J,:)-V(:,1:J-1,:))/da;
        
        %useful quantities
        c_B = max(VbB,10^(-6)).^(-1/ga);
        c_F = max(VbF,10^(-6)).^(-1/ga);
        
        dBB = two_asset_kinked_FOC(VaB,VbB,aaa);
        dFB = two_asset_kinked_FOC(VaB,VbF,aaa);
        dBF = two_asset_kinked_FOC(VaF,VbB,aaa);
        dFF = two_asset_kinked_FOC(VaF,VbF,aaa);
            
        %UPWIND SCHEME
        d_B = (dBF>0).*dBF + (dBB<0).*dBB;

        %state constraints at amin and amax
        d_B(:,1,:) = (dBF(:,1,:)>10^(-12)).*dBF(:,1,:); %make sure d>=0 at amax, don't use VaB(:,1,:)
        d_B(:,J,:) = (dBB(:,J,:)<-10^(-12)).*dBB(:,J,:); %make sure d<=0 at amax, don't use VaF(:,J,:)
        d_B(1,1,:)=max(d_B(1,1,:),0);
        
        %split drift of b and upwind separately
        sc_B = (1-xi)*w*zzz + Rb.*bbb - c_B;
        sd_B = (-d_B - two_asset_kinked_cost(d_B,aaa));
        
        d_F = (dFF>0).*dFF + (dFB<0).*dFB;
        %state constraints at amin and amax
        d_F(:,1,:) = (dFF(:,1,:)>10^(-12)).*dFF(:,1,:); %make sure d>=0 at amin, don't use VaB(:,1,:)
        d_F(:,J,:) = (dFB(:,J,:)<-10^(-12)).*dFB(:,J,:); %make sure d<=0 at amax, don't use VaF(:,J,:)
        
        %split drift of b and upwind separately
        sc_F = (1-xi)*w*zzz + Rb.*bbb - c_F;
        sd_F = (-d_F - two_asset_kinked_cost(d_F,aaa));
        sd_F(I,:,:) = min(sd_F(I,:,:),0);
        
        Ic_B = (sc_B < -10^(-12));
        Ic_F = (sc_F > 10^(-12)).*(1- Ic_B);
        Ic_0 = 1 - Ic_F - Ic_B;
        
        Id_F = (sd_F > 10^(-12));
        Id_B = (sd_B < -10^(-12)).*(1- Id_F);
        Id_B(1,:,:)=0;
        Id_F(I,:,:) = 0; Id_B(I,:,:) = 1; %don't use VbF at bmax so as not to pick up articial state constraint
        Id_0 = 1 - Id_F - Id_B;
        
        c_0 = (1-xi)*w*zzz + Rb.*bbb;
      
        c = c_F.*Ic_F + c_B.*Ic_B + c_0.*Ic_0;
        u = c.^(1-ga)/(1-ga);
        
        %CONSTRUCT MATRIX BB SUMMARING EVOLUTION OF b
        X = -Ic_B.*sc_B/db - Id_B.*sd_B/db;
        Y = (Ic_B.*sc_B - Ic_F.*sc_F)/db + (Id_B.*sd_B - Id_F.*sd_F)/db;
        Z = Ic_F.*sc_F/db + Id_F.*sd_F/db;
        
        for i = 1:Nz
            centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
        end
    
        lowdiag(1:I-1,:) = X(2:I,1,:);
        updiag(2:I,:) = Z(1:I-1,1,:);
        for j = 2:J
            lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];
            updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
        end
        
        for nz=1:Nz
            BBi{nz}=spdiags(centdiag(:,nz),0,I*J,I*J)+spdiags([updiag(:,nz);0],1,I*J,I*J)+spdiags([lowdiag(:,nz);0],-1,I*J,I*J);
        end
        
        BB = [BBi{1}, sparse(I*J,I*J); sparse(I*J,I*J), BBi{2}];
        
        %CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
        dB = Id_B.*dBB + Id_F.*dFB;
        dF = Id_B.*dBF + Id_F.*dFF;
        MB = min(dB,0);
        MF = max(dF,0) + xi*w*zzz + Ra.*aaa;
        MB(:,J,:) = xi*w*zzz(:,J,:) + dB(:,J,:) + Ra(:,J,:).*amax; %this is hopefully negative
        MF(:,J,:) = 0;
        chi = -MB/da;
        yy =  (MB - MF)/da;
        zeta = MF/da;
        
        %MATRIX AAi
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
            
            %Add up the upper, center, and lower diagonal into a sparse matrix
            AAi{nz} = spdiags(AAcentdiag,0,I*J,I*J)+spdiags(AAlowdiag,-I,I*J,I*J)+spdiags(AAupdiag,I,I*J,I*J);
            
        end
        
        AA = [AAi{1}, sparse(I*J,I*J); sparse(I*J,I*J), AAi{2}];
        
        A = AA + BB + Bswitch;
        
        if max(abs(sum(A,2)))>10^(-12)
            disp('Improper Transition Matrix')
            break
        end
        
        
        if max(abs(sum(A,2)))>10^(-9)
           disp('Improper Transition Matrix')
           break
        end
        
        B = (1/Delta + rho)*speye(I*J*Nz) - A;
        
        u_stacked = reshape(u,I*J*Nz,1);
        V_stacked = reshape(V,I*J*Nz,1);
        
        vec = u_stacked + V_stacked/Delta;
        
        V_stacked = B\vec; %SOLVE SYSTEM OF EQUATIONS
            
        V = reshape(V_stacked,I,J,Nz);   
        
        
        Vchange = V - v;
        v = V;
        
        b_dist(:,n)    = max(abs(Vchange(:,:,2)),[],2);
        a_dist(:,n)    = max(max(abs(Vchange),[],3),[],1);
        ab_dist(:,:,n) = max(abs(Vchange),[],3);
       
        dist(n) = max(max(max(abs(Vchange))));
        disp(['Value Function, Iteration ' int2str(n) ', max Vchange = ' num2str(dist(n))]);
        if dist(n)<crit
            disp('Value Function Converged, Iteration = ')
            disp(n)
            break
        end
        
    end
end