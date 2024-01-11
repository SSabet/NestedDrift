function [cpol,dpol] = bdotzero(skipThis,VaF,VaB,a,b,z,I,J,Nz,chi0,chi1,ga,xi,w,ra,rb)

    % Set up objects
    cpol   = zeros(I,J,Nz);
    dpol   = zeros(I,J,Nz);
    dlower = (chi0-1)/chi1;
    
    %%%%% Boundary conditions for amin and amax?
    
    % Looop through the whole state space
    for bi = 1:I
        for aj = 1:J
            for zi = 1:Nz
                
                % Generic test function
                tester = @(x,Va,positive) fone(b(bi),a(ai),z(zi),x,Va,positive,chi0,chi1,ga,xi,w,rb);
                
                % Determine if special region exists
                dupper = -(ra*a(ai) + xi*w*z(zi));
                
                % Policies already identified for this point in the state space i.e. bdot != 0 here
                if skipThis(bi,aj,zi) == 1 
                    continue
                
                % No-arbitrage consistent with positive deposit
                elseif tester(0,VaF(bi,aj,zi),1) < 0 
                    
                    bounds = [0,inf];
                    fun    = @(x) tester(x,VaF(bi,aj,zi),1);
                    d      = fsolve(fun,bounds);
                    
                    dpol(bi,aj,zi) = d;
                    cpol(bi,aj,zi) = rb(bi) * b(bi) + (1-xi) * w * z(zi) - d - two_asset_kinked_cost(d,a(ai));
                
                % There is only one drift direction for d < 0
                elseif dupper < dlower
                    
                    % Inaction region with up-drift
                    if tester(0,VaF(bi,aj,zi),0) <= 0
                        
                        dpol(bi,aj,zi) = 0;
                        cpol(bi,aj,zi) = rb(bi) * b(bi) + (1-xi) * w * z(zi);
                        
                    % Negative d but up-drift region
                    else
                        
                        bounds = [dlower,0];
                        fun    = @(x) tester(x,VaF(bi,aj,zi),0);
                        d      = fsolve(fun,bounds);
                    
                        dpol(bi,aj,zi) = d;
                        cpol(bi,aj,zi) = rb(bi) * b(bi) + (1-xi) * w * z(zi) - d - two_asset_kinked_cost(d,a(ai));
                    end
                
                % Drift can be up or down for d < 0
                else 
                    
                    % Inaction region with up-drift
                    if tester(0,VaF(bi,aj,zi),0) <= 0 
                        
                        dpol(bi,aj,zi) = 0;
                        cpol(bi,aj,zi) = rb(bi) * b(bi) + (1-xi) * w * z(zi);
                        
                    % Negative d but up-drift region
                    elseif tester(dupper,VaF(bi,aj,zi),0) < 0 
                        
                        bounds = [dupper,0];
                        fun    = @(x) tester(x,VaF(bi,aj,zi),0);
                        d      = fsolve(fun,bounds);
                    
                        dpol(bi,aj,zi) = d;
                        cpol(bi,aj,zi) = rb(bi) * b(bi) + (1-xi) * w * z(zi) - d - two_asset_kinked_cost(d,a(ai));
                        
                    % Negative d but down-drift region
                    elseif tester(dupper,VaB(bi,aj,zi),0) > 0
                        
                        bounds = [dlower,dupper];
                        fun    = @(x) tester(x,VaB(bi,aj,zi),0);
                        d      = fsolve(fun,bounds);
                    
                        dpol(bi,aj,zi) = d;
                        cpol(bi,aj,zi) = rb(bi) * b(bi) + (1-xi) * w * z(zi) - d - two_asset_kinked_cost(d,a(ai));
                        
                    % Negative d but no drift
                    else
                       
                        dpol(bi,aj,zi) = dupper;
                        cpol(bi,aj,zi) = rb(bi) * b(bi) + (1-xi) * w * z(zi) - d - two_asset_kinked_cost(d,a(ai));
                        
                    end
                end
            end
        end
    end 
end