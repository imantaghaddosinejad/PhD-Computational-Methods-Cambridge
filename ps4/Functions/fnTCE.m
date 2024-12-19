%% COMPUTE TRANSITION COMPETITIVE EQUILIBRIUM PATH %%%%%%%%%%%%%%%%%%%%%%%%  
%
%   Args:
%       params:
%       opts_fnTCE.Tmax: (scalar) 
%       opts_fnTCE.tfp_path: (vector) TFP path of length Tmax+1
%       opts_fnTCE.aggK_guess: (scalar) initial K guess
%       opts_fnTCE.aggL: (scalar) exog. agregate labour supply 
%       opts_fnTCE.mDist_start: (matrix) t=1 SRCE distribution
%       opts_fnTCE.mDist_end: (matrix) t=Tmax+1 SRCE distribution
%       opts_fnTCE.mVF_start: (matrix) t=1 SRCE VF
%       opts_fnTCE.mVF_end: (matrix) t=Tmax+1 SRCE VF
%       opts_fnTce.pct_check: (scalar) percentage of path to cut and check
%           above for convergence of tail
%       opts_fnTCE.progreport: (boolean) print progress report  
%
%   Returns:
%       SolnPath: (struct) 1x1 struct of all on the TCE path objects
%
function [SolnPath] = fnTCE(params, opts_fnTCE)
    % set maximum time limit 
    T = opts_fnTCE.T;
    
    % create placeholders for TCE objects 
    SolnPath = struct();
    SolnPath.aggK = opts_fnTCE.aggK_guess;
    SolnPath.TFP = opts_fnTCE.tfp_path;
    SolnPath.Kmcc = zeros(T+1, 1);
    SolnPath.aggL = opts_fnTCE.aggL_guess;
    SolnPath.mPolaprime = zeros(params.Na2, params.Nz, T+1);
    SolnPath.mPolc = zeros(params.Na1, params.Nz, T+1);
    SolnPath.mVF = zeros(params.Na1, params.Nz, T+1);
    SolnPath.mCurrDist = zeros(params.Na2, params.Nz, T+1);
   
    % fix endog. capital path at t=1 and t=T+1 (transition to match SRCE) 
    SolnPath.Kmcc(1) = SolnPath.aggK(1);
    SolnPath.Kmcc(T+1) = SolnPath.aggK(T+1);
    
    % set SRCE (t=1) for begining (pre shock) period and end of path (T+1)
    % need V_{T+1} for backward iteration and \Phi_{1} to evolve distribution 
    SolnPath.mCurrDist(:, :, 1) = opts_fnTCE.mDist_start;
    SolnPath.mCurrDist(:, :, T+1) = opts_fnTCE.mDist_end;
    SolnPath.mVF(:, :, 1) = opts_fnTCE.mVF_start;
    SolnPath.mVF(:, :, T+1) = opts_fnTCE.mVF_end;
    
    % initialise figure 
    figure('Name', 'K Path (TCE) for a Temporary TFP Shock')
    title('K Path (TCE) for a 5% Decrease in TFP')
    xlabel('Time')
    ylabel('K')
    grid on;
    
    % TCE loop
    iter = 1;
    err = 10;
    convergedTCE = false;
    MaxIter = 3000;
    errTol = 1e-7;
    tic;
    while err > errTol && iter <= MaxIter
       
        opts_fnVFIBackwardSoln = struct();
        for t = T:(-1):1     
            % solve optimal HH decision for current period 
            opts_fnVFIBackwardSoln.VF = SolnPath.mVF(:, :, t+1);
            opts_fnVFIBackwardSoln.aggK = SolnPath.aggK(t); 
            opts_fnVFIBackwardSoln.aggL = SolnPath.aggL(t); 
            params.tfp = SolnPath.TFP(t);
            [mVF_t, mPolaprime_t, ~] = fnVFIBackwardSoln(params, opts_fnVFIBackwardSoln);
            
            % update
            SolnPath.mVF(:, :, t) = mVF_t;
            SolnPath.mPolaprime(:, :, t) = mPolaprime_t;
        end
        
        % update distribution path 
        opts_fnEvolDist = struct();
        for t = 1:T-1
            opts_fnEvolDist.mCurrDist = SolnPath.mCurrDist(:, :, t);
            opts_fnEvolDist.mPolaprime = SolnPath.mPolaprime(:, :, t);
            SolnPath.mCurrDist(:, :, t+1) = fnEvolDist(params, opts_fnEvolDist);
        end
        
        % update aggregate price path 
        for t = 2:T
            vMarginalDista = sum(SolnPath.mCurrDist(:, :, t), 2);
            SolnPath.Kmcc(t) = vMarginalDista' * params.vGrida2'; % vGrida2 is 1xpNa2 (row) vector
        end
        
        % compute error on transition path 
        err = max(abs(SolnPath.Kmcc - SolnPath.aggK)); % the first and last (T+1) entry in Kmcc and aggK vectors are equal 
        
        % update price path
        SolnPath.aggK = params.wtOld.*SolnPath.aggK + (1-params.wtOld).*SolnPath.Kmcc;
        
        timer = toc;
        if mod(iter, 10) == 0 %%&& opts_fnTCE.progreport
            fprintf('Iteration: %d in %.3f mins. Transition Error %.6f\n', iter, timer/60, err);
            fprintf('----------------------------------------------------\n')
    
            plot(0:T-1, SolnPath.Kmcc(1:T), 'LineStyle', '-', 'LineWidth', 1.2, 'Color', 'b')
            hold on;
            plot(0:T-1, SolnPath.aggK(1:T), 'LineStyle', '--', 'LineWidth', 1.2, 'Color', 'r')
            hold off;
            grid on;
            xlim([0 T-1])
            yline(SolnPath.aggK(T+1),'Color', 'k', 'LineWidth', 0.7, 'LineStyle', '-')
            ylabel('K')
            xlabel('Time')
            title('K Path (TCE)')
            legend('Actual', 'Predicted', 'SS K', 'Location', 'southeast')
            drawnow;
        end
        iter = iter + 1;
    end
    if err <= errTol
        convergedTCE = true;
        fprintf('TCE Converged: %s. After %d iterations.\n', mat2str(convergedTCE), iter)
    end
end