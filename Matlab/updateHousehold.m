function [c,d,sb,sa,vNew,A] = updateHousehold(V,Rb,Ra,d_zerodrift,d_lower,Bswitch,grids,par)

    % Unpack parameters and grids
    cellfun(@(x) assignin('caller', x, par.(x)), fieldnames(par));
    cellfun(@(x) assignin('caller', x, grids.(x)), fieldnames(grids));
    
    aaa = grids.aaa; bbb = grids.bbb; zzz = grids.zzz;
    
    % Preallocate
    VbF = zeros(I,J,Nz);
    VbB = zeros(I,J,Nz);
    VaF = zeros(I,J,Nz);
    VaB = zeros(I,J,Nz);
    c_B = zeros(I,J,Nz);
    c_F = zeros(I,J,Nz);
    
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

    % Idenfity whether the conditional deposit policy creates drift consistent with the conditions used
    I_BB = (d_BB < d_zerodrift);
    I_BF = (d_BF > d_zerodrift);
    I_FF = (d_FF > d_zerodrift);
    I_FB = (d_FB < d_zerodrift);

    % Build deposit policy for backward liquid drift cases, using only consistent deposit policies
    d_B = d_BF.*I_BF + d_BB.*I_BB + d_zerodrift.*(~I_BB .* ~I_BF);
    d_B(1,:,:) = 0; % This case will never be used, just over-writing a nan when updating d
    
    % Build backward liquid drift policy, and an indicator for when it's consitent with itself
    sb_B = (1-xi)*w*zzz + Rb.*bbb - d_B -...
        adjustmentCost(d_B,aaa, chi0, chi1) - c_B;
    
    % At lower b-boundary don't use Vb_B; if Vb_F dosn't work, leave it to the next part that deals with b-drift zero
    I_B = sb_B < 0; 
    I_B(1,:,:) = 0;

    % ...and equivalents for forward liquid drift case
    d_F = d_FF.*I_FF + d_FB.*I_FB + d_zerodrift.*(~I_FB .* ~I_FF);
    d_F(I,:,:) = 0; 
    sb_F = (1-xi)*w*zzz + Rb.*bbb - d_F -...
        adjustmentCost(d_F,aaa, chi0, chi1) - c_F;
    I_F = (sb_F > 0) .* (I_B==0); % Giving precedence to the backward drift if there's a clash
    I_F(I,:,:) = 0;
    
    % Find consumption and deposit policies for the case of zero liquid drift
    I_0 = 1 - I_B - I_F;
    d_0 = bdotzero(I_0,VaF,VaB,b,a,z,Rb,Ra,d_zerodrift,d_lower,par);
    c_0 = (1-xi)*w*zzz + Rb.*bbb - d_0 - adjustmentCost(d_0, aaa, chi0, chi1);
    
    % Build unconditional policies
    c  = c_F.*I_F + c_B.*I_B + c_0.*I_0;
    d  = d_F.*I_F + d_B.*I_B + d_0.*I_0;
    sb = (1-xi)*w*zzz + Rb.*bbb - c - d - adjustmentCost(d, aaa, chi0, chi1);
    sa = Ra.*aaa + xi .* w .* zzz + d;
        
    % if (sum(sum(sum(I_0<0)))>0)
    %     fprintf('There are elements with both I_F & I_B positive. V not concave in b? \n')
    %     [iii, jjj, kkk] =  ind2sub(size(I_0), find(I_0 < 0));
    %     not_ccv_indexes = [find(I_0 < 0), iii, jjj, kkk];    % 
    % end
    % 
    % if (sum(sum(sum(VbB<0)))>0)
    %     fprintf('There are elements with VbB < 0. V not concave in b? \n') 
    %     [iii, jjj, kkk] =  ind2sub(size(VbB), find(VbB < 0));
    %      neg_Vb_indexes = [find(VbB < 0), iii, jjj, kkk];
    % end
    
    % ---------------------------------------------------------------------
    % Value update
    
    % Build transition matrix
    A  = driftMatrixLiquid(sb,bbb,par) + driftMatrixIlliquid(sa,aaa,par) + Bswitch;
    
    % if max(abs(sum(A,2)))>10^(-12)
    % 
    %     disp('Improper Transition Matrix')
    %     [ii, jj, kk] =  ind2sub(size(V), find(abs(sum(A,2))>10^(-12)));
    %     Improper_entries = [find(abs(sum(A,2)) > 10^(-12)), ii, jj, kk];
    % 
    %     break
    % end
    
    % Build update objects
    B         = (1/Delta + rho)*speye(I*J*Nz) - A;
    u_stacked = reshape(U(c,par),I*J*Nz,1);
    V_stacked = reshape(V,I*J*Nz,1);
    vec       = u_stacked + V_stacked/Delta;
    
    % Solve & reshape
    vNew = reshape(B\vec,I,J,Nz);
    
end