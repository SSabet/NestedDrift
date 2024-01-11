function [cpol,dpol] = bdotzero(bdrift_zero,VaF,VaB,a,b,z,Rb,par)
    cellfun(@(x) assignin('caller', x, par.(x)), fieldnames(par));

    % Set up objects
    cpol   = zeros(I,J,Nz);
    dpol   = zeros(I,J,Nz);
    
    
    %%%%% Boundary conditions for amin and amax?
    % for now: looping aj over 2:(J-1), then handling amin and amax seperately
    % after the main aj loop,

    % Loop through the whole state space
    for bi = 1:I
        for aj = 2:(J-1)
            dlower = (chi0-1)*a(aj)/chi1;
            for zk = 1:Nz
                % FOR DEBUGGING
                % if ((bi == 50 || bi == 51) && aj == 2 && zk ==1)
                %     fprintf('cons neg?')
                % end

                % Generic test function
                tester = @(x,Va,positive) fone(b(bi),a(aj),z(zk),x,Va,positive,Rb(bi,aj,zk),par);
                % Determine if special region exists
                dupper = -(ra*a(aj) + xi*w*z(zk));
                
                % Policies already identified for this point in the state space i.e. bdot != 0 here
                if bdrift_zero(bi,aj,zk) == 0
                    continue
                
                % No-arbitrage consistent with positive deposit
                elseif tester(0,VaF(bi,aj,zk),1) < 0 
                    % to search in the positive domain, d cannot be too
                    % large because then consumption would be negative.
                    % let's find a large legitimate d!
                    dmax = max(roots([chi1/(2*a(aj)), chi0, -(Rb(bi,aj,zk)*b(bi) + (1-xi)*w*z(zk))]));

                    bounds = [0,dmax];
                    fun    = @(x) tester(x,VaF(bi,aj,zk),1);
                    d      = fzero(fun,bounds);
                    
                    dpol(bi,aj,zk) = d;
                    cpol(bi,aj,zk) = Rb(bi,aj,zk) * b(bi) + (1-xi) * w * z(zk) - d - two_asset_kinked_cost(d,a(aj), chi0, chi1);
                
                % There is only one drift direction for d < 0
                elseif dupper < dlower
                    
                    % Inaction region with up-drift
                    if tester(0,VaF(bi,aj,zk),0) <= 0
                        
                        dpol(bi,aj,zk) = 0;
                        cpol(bi,aj,zk) = Rb(bi,aj,zk) * b(bi) + (1-xi) * w * z(zk);
                        
                    % Negative d but up-drift region
                    else
                        
                        bounds = [dlower,0];
                        fun    = @(x) tester(x,VaF(bi,aj,zk),0);
                        
                        d = fzero(fun,bounds);
                        
                        dpol(bi,aj,zk) = d;
                        cpol(bi,aj,zk) = Rb(bi,aj,zk) * b(bi) + (1-xi) * w * z(zk) - d - two_asset_kinked_cost(d,a(aj), chi0, chi1);
                    end
                
                % Drift can be up or down for d < 0
                else 
                    
                    % Inaction region with up-drift
                    if tester(0,VaF(bi,aj,zk),0) <= 0 
                        
                        dpol(bi,aj,zk) = 0;
                        cpol(bi,aj,zk) = Rb(bi,aj,zk)* b(bi) + (1-xi) * w * z(zk);
                        
                    % Negative d but up-drift region
                    elseif tester(dupper,VaF(bi,aj,zk),0) < 0 
                        
                        bounds = [dupper,0];
                        fun    = @(x) tester(x,VaF(bi,aj,zk),0);
                        d      = fzero(fun,bounds);
                    
                        dpol(bi,aj,zk) = d;
                        cpol(bi,aj,zk) = Rb(bi,aj,zk) * b(bi) + (1-xi) * w * z(zk) - d - two_asset_kinked_cost(d,a(aj), chi0, chi1);
                        
                    % Negative d but down-drift region
                    elseif tester(dupper,VaB(bi,aj,zk),0) > 0
                        
                        bounds = [dlower,dupper];
                        fun    = @(x) tester(x,VaB(bi,aj,zk),0);
                        if (aj == 1)
                            d = dupper; % at amin we should avoide negative drift
                        else
                            d      = fzero(fun,bounds);
                        end
                    
                        dpol(bi,aj,zk) = d;
                        cpol(bi,aj,zk) = Rb(bi,aj,zk) * b(bi) + (1-xi) * w * z(zk) - d - two_asset_kinked_cost(d,a(aj), chi0, chi1);
                        
                    % Negative d but no drift
                    else
                       
                        dpol(bi,aj,zk) = dupper;
                        cpol(bi,aj,zk) = Rb(bi,aj,zk) * b(bi) + (1-xi) * w * z(zk) - d - two_asset_kinked_cost(dupper,a(aj), chi0, chi1);
                        
                    end
                end
            end
        end
        
        % Handling the amin and amax cases
        % at amin, enforce: d =0
        aj = 1;
        
        for zk = 1:Nz
                dpol(bi,aj,zk) = 0;
                cpol(bi,aj,zk) = Rb(bi,aj,zk) * b(bi) + (1-xi) * w * z(zk);
        end

        % now for amax
        aj = J;
        dlower = (chi0-1)*a(aj)/chi1;
        for zk = 1:Nz
            tester = @(x,Va,positive) fone(b(bi),a(aj),z(zk),x,Va,positive,Rb(bi,aj,zk),par);
            dupper = -(ra*a(aj) + xi*w*z(zk));
            
            if bdrift_zero(bi,aj,zk) == 0
                    continue

            % Negative d but down-drift region: only if dupper > dlower
            % (otherwise drift would be positive: enforce 0)
            elseif (dupper > dlower && tester(dupper,VaB(bi,aj,zk),0) > 0)
                bounds = [dlower,dupper];
                fun    = @(x) tester(x,VaB(bi,aj,zk),0);
                d      = fzero(fun,bounds);
                    
                dpol(bi,aj,zk) = d;
                cpol(bi,aj,zk) = Rb(bi,aj,zk) * b(bi) + (1-xi) * w * z(zk) - d - two_asset_kinked_cost(d,a(aj), chi0, chi1);
            else
                dpol(bi,aj,zk) = dupper;
                cpol(bi,aj,zk) = Rb(bi,aj,zk) * b(bi) + (1-xi) * w * z(zk) - d - two_asset_kinked_cost(dupper,a(aj), chi0, chi1);
            end           
        end
    end
end