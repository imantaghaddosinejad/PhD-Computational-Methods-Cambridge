%% Assignment 4 

clear 
clear all
clc 
cd 'C:\Users\it315\Documents\phd_computational_macro\ps4';
addpath('./Functions')
%addpath('../ReplicationFiles/Aiyagari1994/Matlab')
%addpath('../ReplicationFiles/Aiyagari1994/Matlab/Functions')

%% Ex. 4.a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters 
params.alpha = 0.36;
params.beta = 0.96;
params.delta = 0.08;
params.tfp = 1.0; % initial TFP 
params.mu = 0;
params.rho = 0.9;
params.sigma = 0.2 * sqrt(1 - params.rho^2);
params.Nz = 7;
params.minwealth = 0;
params.initk = 6;
params.progreport = true;

% solve SRCE for t=0 (TFP = 1.0) 
[mPolc_start, mPolaprime_start, mVF_start, mCurrDist_start, ...
    aggK_start, aggL_start, w_start, r_start, ...
    vGrida1, vGrida2, vGridz, mPz, vPiz] = fnSolveSRCE(params);

%% Ex 4.b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.tfp = 1.1; % update TFP (10% increase) 
% solve SRCE for t=T+1 (TFP = 1.1)
[mPolc_end, mPolaprime_end, mVF_end, mCurrDist_end, ...
    aggK_end, aggL_end, w_end, r_end, ...
    ~, ~, ~, ~, ~] = fnSolveSRCE(params);

% plot marginal distributions across different technologies (TFPs)
vMarginalDista_start = sum(mCurrDist_start, 2);
vMarginalDista_end = sum(mCurrDist_end, 2);

plot(vGrida2, vMarginalDista_start, 'LineWidth', 1.8)
hold on;
plot(vGrida2, vMarginalDista_end, 'LineWidth', 1.8)
hold off;
xlabel('Assets')
ylabel('Density')
title('Marginal Distribution of Assets')
grid on;
legend('TFP = 1.0', 'TFP = 1.1')

%% Ex. 4.c, 4.d SOLVING FOR A TRANSITION COMPETITIVE EQUILIBRIUM (TCE) %%%%

% set parameters 
params.tfp = 1.1; % TFP shock at T+1 (10% increase) internalised by agents w/ perfect foresight
params.Na1 = 50;
params.Na2 = 100;
params.mPz = mPz;
params.vGrida1 = vGrida1;
params.vGrida2 = vGrida2;
params.vGridz = vGridz;
params.wtOld = 0.9000;

T = 301; % set transition time limit 

% set key transition objects to begin loop
opts_fnTCE = struct();
opts_fnTCE.T = T; % total number of periods on transition path (set t=1 to initial SRCE). T+1 is terminal SRCE (post shock)
opts_fnTCE.aggL_guess = zeros(T+1, 1);
opts_fnTCE.aggL_guess(1:end) = aggL_end; % exg. labour supply 
opts_fnTCE.mDist_start = mCurrDist_start; % initial (t=1) SRCE distribution over (a,z) 
opts_fnTCE.mDist_end = mCurrDist_end; % terminal (t=T+1) SRCE distribution over (a,z)
opts_fnTCE.mVF_start = mVF_start; % initial (t=1) SRCE VF 
opts_fnTCE.mVF_end = mVF_end; % terminal (t=T+1) SRCE VF
opts_fnTCE.progreport = true; % report progress and plot K path

% set K path - initial guess
opts_fnTCE.aggK_guess = zeros(T+1, 1);
opts_fnTCE.aggK_guess(1) = aggK_start;
opts_fnTCE.aggK_guess(2:end) = aggK_end;

% set TFP path - permanent increase in TFP by 10% 
opts_fnTCE.tfp_path = ones(T+1, 1);
opts_fnTCE. tfp_path(2:T+1) = 1.1;

% compute TCE for permanent shock
SolnPath1 = fnTCE(params, opts_fnTCE);

%% Ex. 4.f SOLVING FOR TRANSITIONAL DYNAMICS UNDER PARTIAL EQUILIBRIUM %%%%

% set parameters 
params.tfp = 1.1; % TFP shock at T+1 (10% increase) internalised by agents w/ perfect foresight

% initialise path for optimal solution (VF, PF and eqlb. prices)
T = 2001; % add one extra period to let t=1 be initial SRCE state (i.e. t=0)
SolnPath.aggK = zeros(T+1, 1); % initial K
SolnPath.Kmcc = SolnPath.aggK; % endog. K computed from market clearing conditions
%SolnPath.aggKnew = SolnPath.aggK; % updated path of capital
SolnPath.aggL = zeros(T+1, 1);
SolnPath.mPolaprime = zeros(params.Na2, params.Nz, T+1);
SolnPath.mPolc = zeros(params.Na1, params.Nz, T+1);
SolnPath.mVF = zeros(params.Na1, params.Nz, T+1);
SolnPath.mCurrDist = zeros(params.Na2, params.Nz, T+1);

% guessed path for aggregate allocation(s) 
for t = 2:T+1
    SolnPath.aggK(t) = aggK_end;
    SolnPath.aggL(t) = aggL_end;
end
SolnPath.aggK(1) = aggK_start; % start SRCE capital allocation
SolnPath.aggL(1) = aggL_start; % start SRCE labour supply allocaiton

% set terminal SRCE equilibrium objects 
SolnPath.mVF(:, :, T+1) = mVF_end; % new SRCE (post TFP shock)
SolnPath.mPolaprime(:, :, T+1) = mPolaprime_end; % new savings policy rule (post TFP shock)
SolnPath.mPolc(:, :, T+1) = mPolc_end; % new consumption policy rule (post TFP shock)
SolnPath.mCurrDist(:, :, T+1) = mCurrDist_end; % new SRCE distribution over (a,z) (post TFP shock)

% initialise joint-distribution (initial SRCE, at t=0) 
SolnPath.mCurrDist(:, :, 1) = mCurrDist_start; % initial SRCE distribution over (a,z) (before TFP shock)

% set the start and end of the endog. K path to match SRCE 
SolnPath.Kmcc(1) = SolnPath.aggK(1);
SolnPath.Kmcc(T+1) = SolnPath.aggK(T+1);

% solve optimal HH problem using backward iteration
figure('Name', 'K Path (Transitional Partial Equilibrium)')
title('K Path (Transitional Partial Equilibrium)')
xlabel('Time')
ylabel('K')
grid on;

err = 10; 
errTol = 1e-6;
iter = 1;
Maxiter = 3000;
%while err > errTol && iter <= Maxiter
while iter < 2
    mVF_t = zeros(params.Na1, params.Nz);
    mPolaprime_t = zeros(params.Na2, params.Nz);
    mPolc_t = zeros(params.Na1, params.Nz);
    
    opts_fnVFIBackwardSoln = struct();
    for t = T:(-1):1 
    
        % solve optimal HH decision for current period 
        opts_fnVFIBackwardSoln.VF = SolnPath.mVF(:, :, t+1); % next period VF (V_{t+1})
        opts_fnVFIBackwardSoln.aggK = SolnPath.aggK(t); % current capital allocation K(t)
        opts_fnVFIBackwardSoln.aggL = SolnPath.aggL(t); % current labour supply L(t)
        [mVF_t, mPolaprime_t, mPolc_t] = fnVFIBackwardSoln(params, opts_fnVFIBackwardSoln);
        
        % update optimal HH policy rules and value function for current period
        SolnPath.mVF(:, :, t) = mVF_t;
        SolnPath.mPolaprime(:, :, t) = mPolaprime_t;
        SolnPath.mPolc(:, :, t) = mPolc_t;
    end
    
    % update distribution path 
    opts_fnEvolDist = struct();
    for t = 1:T-1
        opts_fnEvolDist.mCurrDist = SolnPath.mCurrDist(:, :, t);
        opts_fnEvolDist.mPolaprime = SolnPath.mPolaprime(:, :, t);
        SolnPath.mCurrDist(:, :, t+1) = fnEvolDist(params, opts_fnEvolDist); % update next period dist.
    end
    
    % update aggregate price path 
    for t = 2:T
        vMarginalDista = sum(SolnPath.mCurrDist(:, :, t), 2);
        SolnPath.Kmcc(t) = vMarginalDista' * params.vGrida2'; % vGrida2 is 1xpNa2 (row) vector
    end
    
    % compute error on price path 
    err = max(abs(SolnPath.Kmcc - SolnPath.aggK)); % the first and last (T+1) entry in Kmcc and aggK vectors are equal 

    % update price path
    %SolnPath.aggK = params.wtOld.*SolnPath.aggK + (1-params.wtOld).*SolnPath.Kmcc;
    
    if mod(iter, 1) == 0
        fprintf('Iteration: %d. Error %.9f\n', iter, err);
        fprintf('----------------------------------------------------\n')

        plot(0:T, SolnPath.Kmcc(1:T+1), 'LineStyle', '-', 'LineWidth', 1.2, 'Color', 'b')
        grid on;
        xlim([0 T])
        yline(SolnPath.aggK(T+1),'Color', 'k', 'LineWidth', 0.7, 'LineStyle', '-')
        ylabel('K')
        xlabel('Time')
        title('K Path to a 10% Increase in TFP (Partial Equilibrium)')
        legend('Actual', 'Predicted', 'New SS K', 'Location', 'southeast')
        drawnow;
    end
    iter = iter + 1;
end

%% Ex. 4.h TRANSITIONAL DYNAMICS FOR A TEMPORARY TFP SHOCK %%%%%%%%%%%%%%%%

% set key transition objects to begin loop
T = 201;

% set TCE loop objects
opts_fnTCE = struct();
opts_fnTCE.T = T; 
opts_fnTCE.aggK_guess = aggK_start;
opts_fnTCE.aggL_guess = zeros(T+1, 1);
opts_fnTCE.aggL_guess(1:end) = aggL_end; % exg. labour supply 
opts_fnTCE.mDist_start = mCurrDist_start;
opts_fnTCE.mDist_end = mCurrDist_start;
opts_fnTCE.mVF_start = mVF_start;
opts_fnTCE.mVF_end = mVF_start;
opts_fnTCE.progreport = true;

% set K path - initial guess
opts_fnTCE.aggK_guess = zeros(T+1, 1);
opts_fnTCE.aggK_guess(1:end) = aggK_start;

% set TFP path
opts_fnTCE.tfp_path = zeros(T+1, 1);
opts_fnTCE.tfp_path(1) = 1.0;
opts_fnTCE.tfp_path(2) = 0.95; % temporary (MIT) shock
opts_fnTCE.tfp_path(3:end) = 1.0;

% compute TCE for temporary shock
SolnPath2 = fnTCE(params, opts_fnTCE);

%% Ex 4.i TEMPORARY TFP SHOCK W/ DECAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 301; % set time horizon for Transition path

% set TCE loop objects 
opts_fnTCE = struct();
opts_fnTCE.T = T; 
opts_fnTCE.aggL_guess = zeros(T+1, 1);
opts_fnTCE.aggL_guess(1:end) = aggL_end; 
opts_fnTCE.mDist_start = mCurrDist_start;
opts_fnTCE.mDist_end = mCurrDist_start;
opts_fnTCE.mVF_start = mVF_start;
opts_fnTCE.mVF_end = mVF_start;
opts_fnTCE.progreport = true;

% set K path - initial guess
opts_fnTCE.aggK_guess = zeros(T+1, 1);
opts_fnTCE.aggK_guess(1:end) = aggK_start;

% set TFP Path 
params.rho_tfp = 0.90;
opts_fnTCE.tfp_path = zeros(T+1, 1);
opts_fnTCE.tfp_path(1) = 1.0; % initial SRCE TFP (t=1) 
opts_fnTCE.tfp_path(2) = 0.95; % 5% decrease in TFP 
for t = 3:T+1
    opts_fnTCE.tfp_path(t) = opts_fnTCE.tfp_path(t-1)^params.rho_tfp;
end

% compute TCE for temporary shock
SolnPath3 = fnTCE(params, opts_fnTCE);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SolnPath1.aggKdiff = SolnPath1.aggK - SolnPath1.aggK(T+1);

figure;
plot(SolnPath1.aggKdiff, 'LineWidth', 1)
grid on;
yline(0, 'LineStyle', '--', 'Color', 'k')
xlim([2 90])
ylim([-1 0.2])
xticks(0:10:120)