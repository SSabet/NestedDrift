function sol = solveHousehold(guess,grids,par)

    % Check if there's a starting guess, and generate one if not
    if isnan(guess)
        cGuess = driftLiquid(0,0,grids.bbb,grids.aaa,grids.zzz,par) + ...
            driftilliquid(0,grids.aaa,grids.zzz,par);
        v0 = U(cGuess,par)/par.rho; 
    else
        v0 = guess;
    end
    [vNew,v0_alt] = deal(v0);

    % Solve HJB with tentative guess algorithm
    dist = ones(par.maxit,1); slowCount = 0; fastCount = 1; deltaBase = 0.5;
    for n = 1:par.maxit

        % Relabel the current value guess
        vOld = vNew;
        
        % Tentative guess: attempt an update, and deal with failure if it 
        %   happens by starting again with much sloer update speed
        try

            % Update policies & value with nested-drift algorithm
            [c,d,sb,sa,vNew,A] = updateHousehold(vOld,grids,par);

            % Check convergence & iterate
            dist(n) = max(abs(vOld(:) - vNew(:)));
            disp(['Value Function Iteration ' int2str(n) ', max Vchange = ' num2str(dist(n))]);
            if dist(n) < par.crit
                disp('Value Function Converged, Iteration = ')
                disp(n)
                break
            end 

        % If attempt fails, reset the guess to the last safe starting point
        % and reduce the updating speed to a bare minimum
        catch

            fprintf('\nUpdate failed: looking for better launch point\n')
            vNew      = v0_alt;    % Reset the starting point
            par.Delta = deltaBase; % Slow things down
            slowCount = 1;         % Updating will be slow for 20 iterations

        end
    
        % Tentative guess: control the updating speed
        if slowCount > 0
            slowCount = slowCount + 1;
            if slowCount > 20
                fprintf('\nGoing fast again\n')
                v0_alt    = vNew;             % Set the last slow update as the new safe starting point
                slowCount = 0; fastCount = 1; % Let the update speed start climbing again
            end
        else
            par.Delta = min(1e8,deltaBase*exp(fastCount));
            fastCount = fastCount + 1;
        end
    end

    % Bundle the results into a struct
    sol.c = c; sol.d = d; sol.sb = sb; sol.sa = sa; sol.V = vNew; sol.A = A;

    % Solve the stationary measure
    sol.g  = stationaryDistribution(sol.A, grids.a, grids.b, par);

    % *********************************************************************
    % For interest - check the measure is indeed stationary by rolling the
    % Markov process over a bunch of times
    % g_longrun = (speye(par.Nz*par.I*par.J) - sol.A'*(1000))\sol.g(:); max(abs(g_longrun - sol.g(:)))
    % *********************************************************************

end