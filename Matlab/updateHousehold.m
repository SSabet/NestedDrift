function [c,d,sb,sa,vNew,A] = updateHousehold(V,grids,par)

    % Unpack parameters and grids
    cellfun(@(x) assignin('caller', x, par.(x)), fieldnames(par));
    cellfun(@(x) assignin('caller', x, grids.(x)), fieldnames(grids));
    
    aaa = grids.aaa; bbb = grids.bbb; zzz = grids.zzz;
    
    % Define important policies
    d_zerodrift = - driftilliquid(0,aaa,zzz,par);

    % Preallocate
    [VbF, VbB, VaF, VaB,c_B,c_F] = deal(zeros(I,J,Nz));
    
    % ---------------------------------------------------------------------
    % Policy update
    
    % Prepare value derivatives for use in FOC calculations
    VbF(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))./(bbb(2:I,:,:)-bbb(1:I-1,:,:)); % forward difference wrt b
    VbB(2:I,:,:)   = VbF(1:I-1,:,:); % backward difference wrt b
    
    VaF(:,1:J-1,:) = (V(:,2:J,:)-V(:,1:J-1,:))./(aaa(:,2:J,:)-aaa(:,1:J-1,:)); % forward difference wrt a
    VaB(:,2:J,:)   = VaF(:,1:J-1,:); % backward difference wrt a
        
    % Build consumption and deposit policies conditional on assumed 
    %   liquid and illiquid drift directions 
    c_B(2:I,:,:)   = MU_inv(VbB(2:I,:,:),par);
    c_F(1:I-1,:,:) = c_B(2:I,:,:);

    d_BB = FOC_dPolicy(VaB, VbB, a, par.chi0, par.chi1);
    d_BF = FOC_dPolicy(VaF, VbB, a, par.chi0, par.chi1);
    d_FB = FOC_dPolicy(VaB, VbF, a, par.chi0, par.chi1);
    d_FF = FOC_dPolicy(VaF, VbF, a, par.chi0, par.chi1);

    % Impose the zero drift deposit in cases where the illiquid asset will overrun its state bounds otherwise
    d_FF(:,J,:) = d_zerodrift(:,J,:);
    d_BF(:,J,:) = d_zerodrift(:,J,:);
    d_FB(:,1,:) = d_zerodrift(:,1,:);
    d_BB(:,1,:) = d_zerodrift(:,1,:); 

    % Idenfity whether the conditional deposit policy creates illiquid drift consistent with the conditions used
    I_BB = (d_BB < d_zerodrift);
    I_BF = (d_BF > d_zerodrift);
    I_FF = (d_FF > d_zerodrift);
    I_FB = (d_FB < d_zerodrift);

    % Build deposit policy for backward liquid drift cases, using only consistent deposit policies
    d_B  = d_BF.*I_BF + d_BB.*I_BB + d_zerodrift.*(~I_BB .* ~I_BF);
    d_B(1,:,:) = 0; % This case will never be used, just over-writing a nan when updating d
    
    % Build backward liquid drift policy, and an indicator for when it's consitent with itself
    sb_B = driftLiquid(c_B,d_B,bbb,aaa,zzz,par);
    
    % At lower b-boundary don't use Vb_B; if Vb_F dosn't work, leave it to the next part that deals with b-drift zero
    I_B  = (sb_B < 0); I_B(1,:,:) = 0;

    % ...and equivalents for the forward liquid drift case
    d_F  = d_FF.*I_FF + d_FB.*I_FB + d_zerodrift.*(~I_FB .* ~I_FF); d_F(I,:,:) = 0;
    sb_F = driftLiquid(c_F,d_F,bbb,aaa,zzz,par);
    I_F  = (sb_F > 0) .* (I_B==0); % Giving precedence to the backward drift if there's a clash
    I_F(I,:,:) = 0;
    
    % ...and for the zero liquid drift case
    I_0  = 1 - I_B - I_F;
    d_0  = bdotzero(I_0,VaF,VaB,grids,par);
    c_0  = driftLiquid(0,d_0,bbb,aaa,zzz,par);

    % Combine to build unconditional policies
    c    = c_F.*I_F + c_B.*I_B + c_0.*I_0;
    d    = d_F.*I_F + d_B.*I_B + d_0.*I_0;
    sb   = driftLiquid(c,d,bbb,aaa,zzz,par);
    sa   = driftilliquid(d,aaa,zzz,par);

    % *********************************************************************
    % For interest: show that the zerodrift policy leads to MU(c) between 
    % the forward and backward derivatives (proof in appendix to the paper)
    % mu_check = ((MU(c_0,par) >= VbF) & (MU(c_0,par) <= VbB)).*I_0 + (1-I_0);
    % assert(all(reshape(mu_check(2:I-1,:,:),1,[])), ...
    %     'Error: MU(c_0) is not within the bounds of VbF and VbB at interior points!')
    % *********************************************************************
    
    % ---------------------------------------------------------------------
    % Value update
    
    % Build transition matrix
    A  = driftMatrixLiquid(sb,bbb,par) + driftMatrixIlliquid(sa,aaa,par) + par.Bswitch;
    
    % Build update objects
    B         = (1/Delta + rho)*speye(I*J*Nz) - A;
    u_stacked = reshape(U(c,par),I*J*Nz,1);
    V_stacked = reshape(V,I*J*Nz,1);
    vec       = u_stacked + V_stacked/Delta;
    
    % Solve & reshape
    vNew = reshape(B\vec,I,J,Nz);
    
end
