% Assignment 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;close all;clc;
addpath('./Functions')

%% Part I. 
%%%% Transitional Competitive Equilibrium w/ Permanant TFP Shock %%%%%%%%%%

% parameterize
params.pEta = 6.400;
params.pRiskAversion = 1.000;
params.pBeta = 0.990;
params.pAlpha = 0.330;
params.pDelta = 0.025;
params.pMu = 0.600;
%params.pFrisch = 1.000;

% unpack parameters
pEta=params.pEta;pRiskAversion=params.pRiskAversion;pBeta=params.pBeta;pAlpha=params.pAlpha;pDelta=params.pDelta;pMu=params.pMu;

% compute steady state
ss0.A = 1.000; % initial 
ss0.L = 1/3;
[ss0.L, ss0.K, ss0.C, ss0.Y, ss0.I, ss0.r, ss0.w, pFrisch0] = fnSSCompute(params, ss0,1);

ss1.A = 1.100; % terminal 
ss1.L = 1/3;
[ss1.L, ss1.K, ss1.C, ss1.Y, ss1.I, ss1.r, ss1.w, pFrisch1] = fnSSCompute(params, ss1, 1);
pFrisch=pFrisch1; 

% setup TCE time parameters 
T = 402; % t=1 is initial SS and t=T+1 is terminal SS 
iFuture = [(2:T)';T];
iFuture2 = [(3:T)';T;T];
iPast = [1;(1:T-1)'];

% Shock Path
vA = [ss0.A; ss1.A*ones(T-1,1)]; % shock path (SS(0) is at time t=1, shock kicks in at t=2)
vAFuture = vA(iFuture);

% guess initial solution path 
vK = [ss0.K; ss1.K*ones(T-1,1)]; % anchor t=1 to old steady state 
vC = ss1.C*ones(T,1);
vL = ss1.L*ones(T,1);
vY = ss1.Y*ones(T,1);
vI = ss1.I*ones(T,1);
vr = ss1.r*ones(T,1);
vw = ss1.w*ones(T,1);

% initialize 
vKnew = zeros(T,1);
vCnew = zeros(T,1);
vCfoc = [zeros(T-1,1);ss1.C];

% loop parameters
wtOldK = 0.9000;
wtOldC = 0.9000;
errTol = 1e-5;
MaxIter = 20000;

% TCE loop
iter = 1;
err = 1;
tic;
while err > errTol %&& iter <= MaxIter
    
    % anchor final T+1 value at new SS 
    vK = [vK(1:end-1);ss1.K]; % period T+1 is terminal (new) ss 
    vC = [vC(1:end-1);ss1.C]; % period T+1 is terminal (new) ss 
    vKFuture = vK(iFuture);
    vKFuture2 = vK(iFuture2);
    
    % compute beliefs (RHS of Euler)
    vLFuture = (((1-pAlpha) .* vAFuture .* vKFuture.^pAlpha) ./ (pEta.*vC(iFuture).^pRiskAversion)).^(pFrisch/(1+pFrisch*pAlpha));
    vrFuture = pAlpha .* vAFuture .* vKFuture.^(pAlpha-1) .* vLFuture.^(1-pAlpha) - pDelta;
    
    % given beliefs optimally compute all other intratemporal variables
    vCfoc = (pBeta.*(1+vrFuture - pMu.*((vKFuture2-vKFuture)./vKFuture)./(vKFuture.^2)) ./ vC(iFuture).^pRiskAversion).^(-1/pRiskAversion);
    vL = ((1-pAlpha).*vA.*vK.^pAlpha ./ (pEta.*vCfoc.^pRiskAversion)).^(pFrisch/(1+pFrisch*pAlpha));
    vr = pAlpha .* vA .* vK.^(pAlpha-1) .* vL.^(1-pAlpha) - pDelta;
    vw = (1-pAlpha) .* vA .* vK.^pAlpha .* vL.^(-pAlpha);
    vY = vA .* vK.^pAlpha .* vL.^(1-pAlpha);
    vI = vY - vCfoc - (pMu/2) .* ((vKFuture - vK)./vK).^2;
    
    % Fix period t=1 to become old steady state 
    vKPast = [ss0.K;vK(1:end-1)]; % lagged K values where period t=0 is initial ss 
    vIPast = [ss0.I;vI(1:end-1)]; % lagged I values where period t=0 is initial ss     
    
    % Update K-path (and C-path)
    vKnew = (1-pDelta).*vKPast + vIPast;
    vCnew = vA.*vKnew.^pAlpha.*vL.^(1-pAlpha) - vI - (pMu/2) .* ((vKFuture - vK)./vK).^2;
  
    % compute error 
    err = max(abs(vK - vKnew));
    %err2 = max(abs(vK(1:end-1) - vKnew(1:end-1))); % error w/o tail of K-path
    MSE_C = mean((vC - vCnew).^2);
    MSE_K = mean((vK - vKnew).^2);

    % update K path (TCE)
    vK = wtOldK.*vK + (1-wtOldK).*vKnew;
    vC = wtOldC.*vC + (1-wtOldC).*vCnew;
    
    % print progress 
    timer = toc;
    if mod(iter, 50) == 0
        fprintf('Iteration %d after %.2fs. K Error: %.10f\n', iter, timer, err);    
        fprintf('MSE_C: %.8f. MSE_K: %.8f.\n', MSE_C, MSE_K);
        fprintf('--------------------------------------------------\n');
        
        plot(1:T, vKnew, 'b-', 'LineWidth', 1.1);hold on;
        yline(ss1.K, 'k--', 'Label','New SS');
        plot(1:T, vK, 'r--', 'LineWidth', 0.9); hold off;xlim([1, T]);
        grid on;legend('Actual', 'Guess', '', 'Location','southeast');ylabel('K');drawnow
    end  
    iter = iter+1;
end

% save optimal paths 
SolnPath1.vK=vK;SolnPath1.vC=vC;SolnPath1.vI=vI;SolnPath1.vL=vL;SolnPath1.vY=vY;SolnPath1.vr=vr;SolnPath1.vw=vw;

%%%% Transitional Competitive Equilibrium w/ Temporary TFP Shock %%%%%%%%%%

% Shock Path 
vt = 1:T;
vA = [ss0.A;0.95*ss0.A; zeros(T-2,1)];
for t = 3:T 
    vA(t) = vA(t-1)^0.95; % shock w/ decay 
end
vAFuture = vA(iFuture);

% guess initial solution path 
vK = ss0.K * ones(T,1);
vC = ss0.C*ones(T,1);
vL = ss0.L*ones(T,1);
vY = ss0.Y*ones(T,1);
vI = ss0.I*ones(T,1);
vr = ss0.r*ones(T,1);
vw = ss0.w*ones(T,1);

% initialize 
vKnew = zeros(T,1);
vCnew = zeros(T,1);
vCfoc = zeros(T,1);

% loop parameters
wtOldK = 0.9000;
wtOldC = 0.9000;
errTol = 1e-5;
MaxIter = 20000;

% TCE loop
iter = 1;
err = 1;
tic;
while err > errTol %&& iter <= MaxIter
    
    % anchor final T+1 value at new SS 
    vK = [vK(1:end-1);ss0.K]; % period T+1 is initial SS 
    vC = [vC(1:end-1);ss0.C]; % period T+1 is initial SS 
    vKFuture = vK(iFuture);
    vKFuture2 = vK(iFuture2);

    % compute beliefs (RHS of Euler)
    vLFuture = (((1-pAlpha) .* vAFuture .* vKFuture.^pAlpha) ./ (pEta.*vC(iFuture).^pRiskAversion)).^(pFrisch/(1+pFrisch*pAlpha));
    vrFuture = pAlpha .* vAFuture .* vKFuture.^(pAlpha-1) .* vLFuture.^(1-pAlpha) - pDelta;
    
    % given beliefs optimally compute all other intratemporal variables
    vCfoc = (pBeta.*(1+vrFuture - pMu.*((vKFuture2-vKFuture)./vKFuture)./(vKFuture.^2)) ./ vC(iFuture).^pRiskAversion).^(-1/pRiskAversion);
    vL = ((1-pAlpha).*vA.*vK.^pAlpha ./ (pEta.*vCfoc.^pRiskAversion)).^(pFrisch/(1+pFrisch*pAlpha));
    vr = pAlpha .* vA .* vK.^(pAlpha-1) .* vL.^(1-pAlpha) - pDelta;
    vw = (1-pAlpha) .* vA .* vK.^pAlpha .* vL.^(-pAlpha);
    vY = vA .* vK.^pAlpha .* vL.^(1-pAlpha);
    vI = vY - vCfoc - (pMu/2) .* ((vKFuture - vK)./vK).^2;
    
    % Fix period t=1 to become old steady state 
    vKPast = [ss0.K;vK(1:end-1)]; % lagged K values where period t=0 is initial ss 
    vIPast = [ss0.I;vI(1:end-1)]; % lagged I values where period t=0 is initial ss     
    
    % Update K-path (and C-path)
    vKnew = (1-pDelta).*vKPast + vIPast; 
    vCnew = vA.*vKnew.^pAlpha.*vL.^(1-pAlpha) - vI - (pMu/2) .* ((vKFuture - vK)./vK).^2;
  
    % compute error 
    err = max(abs(vK - vKnew));
    %err2 = max(abs(vK(1:end-1) - vKnew(1:end-1))); % error w/o tail of K-path
    MSE_C = mean((vC - vCnew).^2);
    MSE_K = mean((vK - vKnew).^2);

    % update K path (TCE)
    vK = wtOldK.*vK + (1-wtOldK).*vKnew;
    vC = wtOldC.*vC + (1-wtOldC).*vCnew;
    
    % print progress 
    timer = toc;
    if mod(iter, 100) == 0
        fprintf('Iteration %d after %.2fs. K Error: %.10f\n', iter, timer, err);    
        fprintf('MSE_C: %.8f. MSE_K: %.8f.\n', MSE_C, MSE_K);
        fprintf('--------------------------------------------------\n');
        
        plot(1:T, vKnew, 'b-', 'LineWidth', 1.1);hold on;
        yline(ss0.K, 'k--', 'Label','New SS');
        plot(1:T, vK, 'r--', 'LineWidth', 0.9); hold off;xlim([1, T]);
        grid on;legend('Actual', 'Guess', '', 'Location','southeast');ylabel('K');drawnow
    end  
    iter = iter+1;
end

% save optimal paths 
SolnPath2.vK=vK;SolnPath2.vC=vC;SolnPath2.vI=vI;SolnPath2.vL=vL;SolnPath2.vY=vY;SolnPath2.vr=vr;SolnPath2.vw=vw;

%% Part II. 

% GLOBAL NONLINEAR SOLUTION TO RBC w/ AGG. UNCERTAINTY/SHOCK %%%%%%%%%%%%%%

% TFP parameters 
params.pNA = 7;
params.pRho = 0.95;
params.pSigmaEps = 0.009;
[vAstates, mPA] = fnTauchen(params.pRho, params.pSigmaEps, 0.0, params.pNA);
vGridA = exp(vAstates);

% set time
Tmax = 10000;
BurnT = 500;
T = Tmax + BurnT;
iFuture = [(2:T)';T]; % lead 1 period ahead 
iFuture2 = [(3:T)';T;T]; % lead 2 periods ahead 

% simulate TFP shock 
ivA = fnSimShock(mPA,T,2,4321);
ivAFuture = ivA(iFuture);
vA = vGridA(ivA);
vAFuture = vGridA(ivAFuture);
vAPrRealized = zeros(T,1);
for t = 1:T
    vAPrRealized(t) = mPA(ivA(t), ivAFuture(t));
end

% guess solution path (initial) 
vK = ss0.K + normrnd(0, 0.000001, T, 1); % anchor variable (initial guess all other paths)
vC = ss0.C .* ones(T, 1); % anchor variable (initial guess pins all other paths)
vL = ((1-pAlpha).*vA.*vK.^pAlpha ./ (pEta.*vC.^pRiskAversion)).^(pFrisch / (1 + pFrisch*pAlpha));
vr = pAlpha.*vA.*(vK./vL).^(pAlpha-1) - pDelta;
vw = (1-pAlpha).*vA.*(vK/vL).^pAlpha;
vY = vA.*vK.^pAlpha.*vL.^(1-pAlpha);
vI = vY - vC - (pMu/2).*( (vK(iFuture)-vK)./vK ).^2;

% placeholder variables 
vKnew = zeros(T,1);
vCnew = ss0.C * ones(T,1);

% loop parameters 
errTol = 1e-6;
MaxIter = 20000;
wtOldK = 0.9000;
wtOldC = 0.9000;

% run loop (solve model)
iter = 1;
err = 10;
tic;
while err > errTol && iter <= MaxIter
    
    % ==============================
    % STEP 1. BACKWARD SOLVE FOR BELIEFS 
    % ==============================

    % update future endog. state 
    vKFuture = [vK(2:end);vK(1)]; % anchor off the path value T+1 to initial value (which will always converge to SS-K
    vKFuture2 = [vK(3:end);vK(1);vK(1)];
    
    % compute beliefs over counterfactual states (unrealized)
    vBeliefV1 = 0;
    for iA = 1:params.pNA
        
        % counterfactual exog. shock
        AFuture = vGridA(iA); 
        
        % location (time period) of when counterfactual was realized 
        vCanLoc = find(vA == AFuture);
        vCanLoc(vCanLoc < BurnT) = [];
        vCanLoc(vCanLoc > T-BurnT) = [];

        % K-value candidates for counterfactual, i.e. (K',A') states
        vCan = vK(vCanLoc);
        [vCan, index] = sort(vCan);
        vCanLoc = vCanLoc(index);

        % interpolate K' on K-candidates 
        nLow = sum(repmat(vCan',T,1) < vKFuture, 2);
        nLow(nLow<=1) = 1;
        nLow(nLow>=length(index)) = length(index) - 1;
        nHigh = nLow + 1;
        wtLow = (vK(nHigh) - vKFuture)./(vK(nHigh)-vK(nLow));
        wtLow(wtLow>1) = 1;
        wtLow(wtLow<0) = 0;
        wtHigh = 1 - wtLow;

        % linearly interpolate beliefs (RHS of Euler) given counterfactual
        vLLow = ((1-pAlpha) * AFuture .* vKFuture.^pAlpha ./ (pEta.*vC(vCanLoc(nLow)).^pRiskAversion)).^(pFrisch / (1 + pFrisch*pAlpha));
        vLHigh = ((1-pAlpha) * AFuture .* vKFuture.^pAlpha ./ (pEta.*vC(vCanLoc(nHigh)).^pRiskAversion)).^(pFrisch / (1 + pFrisch*pAlpha));
        vrLow = pAlpha * AFuture .* (vKFuture./vLLow).^(pAlpha-1) - pDelta;
        vrHigh = pAlpha * AFuture .* (vKFuture./vLHigh).^(pAlpha-1) - pDelta;
        
        % cumulatively update beliefs for counterfactuals 
        vBeliefV1 = vBeliefV1 + ... 
            (ivAFuture ~= iA) .* ... % only update belief for counterfactuals (unrealised states tomorrow)
            pBeta .* mPA(ivA, iA) .* ...
            ( wtLow.*(1 + vrLow - (pMu./vKFuture.^2).*((vKFuture2-vKFuture)./vKFuture))./vC(vCanLoc(nLow)).^pRiskAversion + ...
               wtHigh.*(1 + vrHigh - (pMu./vKFuture.^2).*((vKFuture2-vKFuture)./vKFuture))./vC(vCanLoc(nHigh)).^pRiskAversion );
    end
    
    % cumulatively update beliefs for actual realised future states 
    vLFuture = ((1-pAlpha) .* vAFuture .* vKFuture.^pAlpha ./ (pEta.*vC(iFuture).^pRiskAversion)).^(pFrisch / (1 + pFrisch*pAlpha));
    vrFuture = pAlpha .* vAFuture .* (vKFuture./vLFuture).^(pAlpha-1) - pDelta;
    
    vBeliefV1 = vBeliefV1 + ... 
        pBeta .* vAPrRealized .* ... % consider all realised future states 
        ((1 + vrFuture - (pMu./vKFuture.^2).*((vKFuture2-vKFuture)./vKFuture)) ./ vC(iFuture).^pRiskAversion);
    
    % ==============================
    % STEP 2. GIVEN BELIEFS OPTIMALLY SOLVE (BACKWARD) FOR CONTROLS
    % ==============================
    
    vCfoc = (1./vBeliefV1).^(1/pRiskAversion);
    vL = (((1-pAlpha).*vA.*vK.^pAlpha) ./ (pEta.*vCfoc.^pRiskAversion)).^(pFrisch/(1+pFrisch*pAlpha));
    vr = pAlpha.*vA.*(vK./vL).^(pAlpha-1) - pDelta;
    vw = (1-pAlpha).*vA.*(vK./vL).^pAlpha;
    vY = vA.*vK.^pAlpha.*vL.^(1-pAlpha);
    vI = vY - vCfoc - (pMu/2).*((vKFuture-vK)./vK).^2;

    % ==============================
    % STEP 3. SIMULATE FORWARD OPTIMAL NONLINEAR PATH OF CAPITAL
    % ==============================

    vIPast = [ss0.I;vI(1:end-1)]; % lagged I anchored on SS value for t=0 
    vKPast = [ss0.K;vK(1:end-1)]; % lagged K anchored on SS value for t=0
    vKnew = (1-pDelta).*vKPast + vIPast;
    vCnew = vA.*vKnew.^pAlpha.*vL.^(1-pAlpha) - vI - (pMu/2).*((vKnew(iFuture)-vKnew)./vKnew).^2;

    % ==============================
    % COMPUTE MSE FOR UPDATED PATH
    % ==============================

    MSE_C = mean((vC - vCnew).^2);
    MSE_K = mean((vK - vKnew).^2);
    err = mean([vC-vCnew; vK-vKnew].^2);

    % ==============================
    % UPDATE (C,K) PATH
    % ==============================
    
    vK = wtOldK.*vK + (1-wtOldK).*vKnew;
    vC = wtOldC.*vC + (1-wtOldC).*vCnew;

    % ==============================
    % PROGRESS REPORTING 
    % ==============================
    
    timer = toc;
    if mod(iter, 100)==0
        fprintf('Iteration %d. after %.2fs. MSE: %.10f\n', iter, timer, err);
        fprintf('MSE_C: %.6f. MSE_K: %.6f\n', MSE_C, MSE_K);
        fprintf('----------------------------------------\n')
        
        % plots 
        subplot(1,2,1);
        plot(1:T, vKnew(1:T), 'b-', 'LineWidth', 1.1);hold on;
        plot(1:T, vK(1:T), 'r-.', 'LineWidth', .9);
        yline(ss0.K, 'k--', 'LineWidth', 1, 'Label', 'SS K');hold off;
        grid on;xlabel('Time');ylabel('K');xlim([1,T])
        legend('Actual', 'Predicted', '', 'Location', 'northeast')
        
        subplot(1,2,2);
        plot(1:T, vCnew(1:T), 'b-', 'LineWidth', 1.1);hold on;
        plot(1:T, vC(1:T), 'r-.', 'LineWidth', .9);
        yline(ss0.C, 'k--', 'LineWidth', 1, 'Label', 'SS C');hold off;
        grid on;xlabel('Time');ylabel('C');xlim([1,T])
        legend('Actual', 'Predicted', '', 'Location', 'northeast')
        drawnow;pause(0.2);    
    end
    iter = iter + 1;
end



