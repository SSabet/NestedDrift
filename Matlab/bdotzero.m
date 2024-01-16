function dpol = bdotzero(bdrift_zero,VaF,VaB,a,b,z,Rb,Ra, d_upper, d_lower, MU, par)
    
    cellfun(@(x) assignin('caller', x, par.(x)), fieldnames(par));

    % Set up solution object
    dpol = zeros(I,J,Nz);

    % Loop through the whole state space
    for zk = 1:Nz
        for aj = 1:(J-1)
            
            dlower = d_lower(1, aj, zk);

            for bi = 1:I

                % Generic test function & state-specific zero drift policy
                tester = @(x,Va,positive) FOC_bdotzero_resid(b(bi),a(aj),z(zk),x,Va,positive,Rb(bi,aj,zk), MU, par);
                dupper = d_upper(bi,aj,zk);

                % Policies already identified for this point in the state space i.e. bdot != 0 here
                if bdrift_zero(bi,aj,zk) == 0
                    continue

                % No-arbitrage consistent with positive deposit (d>0)
                elseif tester(0,VaF(bi,aj,zk),1) < 0

                    % d cannot be too large because c would need to be
                    % negative to ensure zero liquid drift. Find the max
                    dmax = max(roots([chi1/(2*a(aj)), chi0+1, -(Rb(bi,aj,zk)*b(bi) + (1-xi)*w*z(zk))]))-1e-12;

                    bounds = [0,dmax];
                    fun    = @(x) tester(x,VaF(bi,aj,zk),1);
                    dpol(bi,aj,zk) = fzero(fun,bounds);

                % d <= 0 and illiquid drift always up
                elseif dupper < dlower

                    % d = 0 with illiquid up-drift
                    if tester(0,VaF(bi,aj,zk),0) <= 0

                        dpol(bi,aj,zk) = 0;
                        
                    % d < 0 with illiquid up-drift
                    elseif tester(dlower,VaF(bi,aj,zk),0) <= 0
                        bounds = [dlower,0];
                        fun    = @(x) tester(x,VaF(bi,aj,zk),0);
                        dpol(bi,aj,zk) = fzero(fun,bounds);
                    else
                        dpol(bi,aj,zk) = dlower; % This should never happen (marginal destruction of value, will only seem right if VaB < 0)
                    end

                % d <= 0 and illiquid drift can be up, down or stationary
                else 

                    % d = 0 with illiquid up-drift (note with d=0 down-drift 
                    %   is impossible, so no other cases need be considered)
                    if tester(0,VaF(bi,aj,zk),0) <= 0 
                        dpol(bi,aj,zk) = 0;

                    % d < 0 with illiquid up-drift
                    elseif tester(dupper,VaF(bi,aj,zk),0) <= 0 

                        bounds = [dupper,0];
                        fun    = @(x) tester(x,VaF(bi,aj,zk),0);
                        dpol(bi,aj,zk) = fzero(fun,bounds);

                    % d < 0 with illiquid down-drift
                    elseif tester(dupper,VaB(bi,aj,zk),0) > 0

                        if (aj == 1)
                            dpol(bi,aj,zk) = dupper; % Impose state-constraint (no negative drift) at lower illiquid bound
                        elseif tester(dlower,VaB(bi,aj,zk),0) <= 0
                            bounds = [dlower,dupper];
                            fun    = @(x) tester(x,VaB(bi,aj,zk),0);
                            dpol(bi,aj,zk) = fzero(fun,bounds);
                        else
                            dpol(bi,aj,zk) = dlower; % This should never happen (marginal destruction of value, will only seem right if VaB < 0)
                        end

                    % d < 0 with illiquid zero-drift
                    else
                        dpol(bi,aj,zk) = dupper;
                    end
                end
            end
        end
       
        % Deal with the top end of the illqiuid grid separately
        aj     = J;
        dlower = (chi0-1)*a(aj)/chi1;
        
        for bi = 1:I
            
            tester = @(x,Va,positive) FOC_bdotzero_resid(b(bi),a(aj),z(zk),x,Va,positive,Rb(bi,aj,zk), MU,par);
            dupper = -(Ra(bi,aj,zk)*a(aj) + xi*w*z(zk));

            % Only consider points in the state space without an already
            % defined policy
            if bdrift_zero(bi,aj,zk) == 0
                continue

            % d < 0 and illiquid down-drift: only if dupper > dlower
            % (otherwise drift would be positive: enforce 0)
            elseif (dupper > dlower && tester(dupper,VaB(bi,aj,zk),0) > 0)
                
                if tester(dlower,VaB(bi,aj,zk),0) <= 0
                    bounds = [dlower,dupper];
                    fun    = @(x) FOC_bdotzero_resid(b(bi),a(aj),z(zk),x,VaB(bi,aj,zk),0,Rb(bi,aj,zk), MU,par);
                    dpol(bi,aj,zk) = fzero(fun,bounds);
                else
                    dpol(bi,aj,zk) = dlower; % This should never happen (marginal destruction of value, will only seem right if VaB < 0)
                end
              
            % d < 0 and zero illiquid drift
            else
                
                dpol(bi,aj,zk) = dupper; % This forces a sub-optimal condition to respect the boundary
                
            end           
        end
    end
end
