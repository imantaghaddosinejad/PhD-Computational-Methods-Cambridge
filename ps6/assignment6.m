% Assignment 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('./Functions')
addpath('./Figures')

% parameterize 
params.pBeta = 0.96;
params.pRiskAversion = 1.000;
params.pFrisch = 1.000;
params.pAlpha = 0.330;
params.pDelta = 0.100;
params.pPhi = 0.975;
%params.pEta = 6.400;
params.pRho = 0.900;
params.pSigmaEps = 0.013;
params.pNA = 7;

% unpack parameters 
pBeta = params.pBeta; pRiskAversion = params.pRiskAversion; pFrisch = params.pFrisch; pAlpha = params.pAlpha; pDelta = params.pDelta; pPhi = params.pPhi;

%% compute SS 
ss.A = 1.000;
ss.L = 1/3; % calibrate pEta to match 1/3 SS hours 
ss.r = 1/pBeta - 1;

pEta = (ss.A * (1-pAlpha) * ((ss.r+pDelta)/(pAlpha*ss.A))^((pAlpha-pRiskAversion)/(pAlpha-1)) * ((ss.r+pDelta)/pAlpha - pDelta)^(-pRiskAversion)) /...
    ss.L^(1/pFrisch + pRiskAversion);
ss.K = ((ss.r+pDelta) / (pAlpha*ss.A))^(1/(pAlpha-1)) * ss.L;
ss.C = ((ss.r+pDelta)/pAlpha - pDelta) * ((ss.r+pDelta)/(pAlpha*ss.A))^(1/(pAlpha-1)) * ss.L;
ss.Y = ss.A * ((ss.r+pDelta)/(pAlpha*ss.A))^(pAlpha/(pAlpha-1)) * ss.L;
ss.w = (1-pAlpha) * ss.A * ((ss.r+pDelta)/(pAlpha*ss.A))^(pAlpha/(pAlpha-1));
ss.I = pDelta*ss.K;

%% Global Nonlinear Solution in Sequence Space - Repeated Transition Method

% time parameters 
Tmax = 10000;
BurnT = 500;
T = Tmax + BurnT;
iFuture = [(2:T)';T];

% simulate shock (TFP) path
[vGridA, mPA] = fnTauchen(params.pRho,params.pSigmaEps,0.0,params.pNA);
vGridA = exp(vGridA);
ivA = fnSimShock(mPA, T, 1, 200);
vA = vGridA(ivA);
ivAFuture = ivA(iFuture);
vAFuture = vA(iFuture);
vPARealized = zeros(T,1);
for t = 1:T
    vPARealized(t) = mPA(ivA(t), ivAFuture(t));
end

% guess solution path 
vK = ss.K + normrnd(0, 0.0000001, T, 1);
vC = ss.C .* ones(T,1);
vLambda = zeros(T,1); % initital occ binding guess is that is never binds

% placeholder 
vKnew = zeros(T,1);
vCnew = zeros(T,1);
vLambdanew = zeros(T,1);

% optimal allocaitons path (initial)
vL = ss.L .* ones(T,1);
vr = ss.r .* ones(T,1);
vw = ss.w .* ones(T,1);
vY = ss.Y .* ones(T,1);
vI = ss.I .* ones(T,1);

% loop parameters 
wtOldK      = 0.9000;
wtOldC      = 0.9000;
wtOldLambda = 0.9000;
errTol = 1e-7;
MaxIter = 30000;

% run loop (solve model) 
iter = 1;
err1 = 10;
err2 = 10;
tic;
while err2 > errTol && iter <= MaxIter

    % ============================================================
    % STEP 1. BACKWARD SOLVE FOR BELIEFS 
    % ============================================================

    vKFuture = [vK(2:end);vK(1)];
    vBelief = 0;
    for iA = 1:params.pNA        

        % location (time period) of when counterfactual was realized 
        AFuture = vGridA(iA);
        vCanLoc = find(ivA == iA);
        vCanLoc(vCanLoc > T-BurnT) = [];
        vCanLoc(vCanLoc < BurnT) = [];

        % K-value candidates for counterfactual, i.e. (K',A') states
        vCan = vK(vCanLoc);
        [vCan, index] = sort(vCan);
        vCanLoc = vCanLoc(index);

        % interpolate K' on K-candidates 
        nLow = sum(repmat(vCan', length(vK), 1) < vKFuture, 2);
        nLow(nLow <= 1) = 1;
        nLow(nLow >= length(index)) = length(index) - 1;
        nHigh = nLow + 1;
        wtLow = (vK(nHigh) - vKFuture) ./ (vK(nHigh) - vK(nLow));
        wtLow(wtLow > 1) = 1;
        wtLow(wtLow < 0) = 0;
        wtHigh = 1 - wtLow;

        % linearly interpolate beliefs (RHS of Euler) given counterfactual       
        vLambdaPrime = wtLow.*vLambda(vCanLoc(nLow)) + wtHigh.*vLambda(vCanLoc(nHigh));
        vCLow = vC(vCanLoc(nLow));
        vCHigh = vC(vCanLoc(nHigh));
        vLLow = ((1-pAlpha) * AFuture .* vKFuture.^pAlpha ./ (pEta .* vCLow.^pRiskAversion)).^(pFrisch / (1+pFrisch*pAlpha));
        vLHigh = ((1-pAlpha) * AFuture .* vKFuture.^pAlpha ./ (pEta .* vCHigh.^pRiskAversion)).^(pFrisch / (1+pFrisch*pAlpha));        
        vrLow = pAlpha .* AFuture .* (vKFuture./vLLow).^(pAlpha-1) - pDelta;
        vrHigh = pAlpha .* AFuture .* (vKFuture./vLHigh).^(pAlpha-1) - pDelta;

        % cumulatively update beliefs for counterfactuals
        vBelief = vBelief + ...
            (ivAFuture ~= iA) .* ...
            pBeta .* mPA(ivA, iA) .* ...
            ( wtLow.*(1+vrLow)./vCLow.^pRiskAversion + wtHigh.*(1+vrHigh)./vCHigh.^pRiskAversion - vLambdaPrime.*(1-pDelta));                    

    end
    % cumulatively update beliefs for actual realised future states 
    vLFuture = ((1-pAlpha) * vAFuture .* vKFuture.^pAlpha ./ (pEta .* vC(iFuture).^pRiskAversion)).^(pFrisch / (1+pFrisch*pAlpha));
    vrFuture = pAlpha .* vAFuture .* (vKFuture./vLFuture).^(pAlpha-1) - pDelta;
    
    vBelief = vBelief + ... 
        pBeta .* vPARealized .* ... % consider all realised future states 
        ((1+vrFuture)./vC(iFuture).^pRiskAversion - vLambda(iFuture).*(1-pDelta)); 

    % ============================================================
    % STEP 2. GIVEN BELIEFS OPTIMALLY SOLVE (BACKWARD) FOR ALLOCATIONS
    % ============================================================    
    
    vLambdanew = 1./vC.^pRiskAversion - vBelief;
    vCfoc = (vBelief + vLambda).^(-1/pRiskAversion);
    vL = ((1-pAlpha) * vA .* vK.^pAlpha ./ (pEta .* vCfoc.^pRiskAversion)).^(pFrisch / (1+pFrisch*pAlpha));
    vr = pAlpha .* vA .* (vK./vL).^(pAlpha-1) - pDelta;
    vw = (1-pAlpha) .* vA .* (vK./vL).^pAlpha;
    vY = vA .* vK.^pAlpha .* vL.^(1-pAlpha);
    vI = vY - vCfoc;

    % impose friction on investment (lower bound on investing)
    vI(vI <= pPhi*ss.I) = pPhi*ss.I;
    vLambdanew(vI > pPhi*ss.I) = 0;

    % ============================================================
    % STEP 3. SIMULATE FORWARD
    % ============================================================

    % simulate aggregate capital path (endog. state) and consumption path
    vIPast = [ss.I;vI(1:end-1)]; % lagged I anchored on SS value for t=0 
    vKPast = [ss.K;vK(1:end-1)]; % lagged K anchored on SS value for t=0
    vKnew = (1-pDelta).*vKPast + vIPast;
    vCnew = vA .* vKnew.^pAlpha .* vL.^(1-pAlpha) - vI;

    % ============================================================
    % STEP 4. COMPUTE ERROR FOR UPDATED PATH
    % ============================================================
    
    % pointwise error 
    err1 = max([abs(vK - vKnew); abs(vC - vCnew); abs(vLambda - vLambdanew)]);

    % MSE 
    MSE_C = mean((vC - vCnew).^2);
    MSE_K = mean((vK - vKnew).^2);
    MSE_Lambda = mean((vLambda - vLambdanew).^2);
    err2 =  mean(([vC-vCnew; vK-vKnew; vLambda-vLambdanew]).^2);

    % ============================================================
    % STEP 5. UPDATING
    % ============================================================

    vK      = wtOldK.*vK +              (1-wtOldK).*vKnew;
    vC      = wtOldC.*vC +              (1-wtOldC).*vCnew;
    vLambda = wtOldLambda.*vLambda +    (1-wtOldLambda).*vLambdanew;

    % ============================================================
    % PROGRESS REPORTING 
    % ============================================================
    
    timer = toc;
    if mod(iter, 100)==0
        fprintf('Iteration %d. after %.2fs. MSE: %.12f\n', iter, timer, err2);
        fprintf('MSE_C: %.6f. MSE_K: %.6f. MSE_Lambda: %.6f\n', MSE_C, MSE_K, MSE_Lambda);
        fprintf('Pointwise Error: %.6f\n', err1);
        fprintf('----------------------------------------\n')
        
        % plots 
        subplot(2,2,1);
        plot(1:T, vKnew(1:T), 'b-', 'LineWidth', .8);hold on;
        plot(1:T, vK(1:T), 'r-.', 'LineWidth', .8);
        yline(ss.K, 'k--', 'LineWidth', 1, 'Label', 'SS K');hold off;
        grid on;xlabel('Time');ylabel('K');xlim([1,T])
        legend('Actual', 'Predicted', '', 'Location', 'northeast')
        
        subplot(2,2,2);
        plot(1:T, vCnew(1:T), 'b-', 'LineWidth', .8);hold on;
        plot(1:T, vC(1:T), 'r-.', 'LineWidth', .8);
        yline(ss.C, 'k--', 'LineWidth', 1, 'Label', 'SS C');hold off;
        grid on;xlabel('Time');ylabel('C');xlim([1,T])
        legend('Actual', 'Predicted', '', 'Location', 'northeast')

        subplot(2,2,3);
        plot(1:T, vLambdanew(1:T), 'b-', 'LineWidth', .8);hold on;
        plot(1:T, vLambda(1:T), 'r-.', 'LineWidth', .8);
        yline(0, 'k--', 'LineWidth', 1, 'Label', 'SS \lambda');hold off;
        grid on;xlabel('Time');ylabel('\lambda');xlim([1,T])
        legend('Actual', 'Predicted', '', 'Location', 'northeast')        
        drawnow;pause(0.2);    

    end
    iter = iter + 1;

    % control how smoothly prices update conditional on iteration number 
    % as iteration blows up, slow down (increase smoothness) price updating
    % this avoids jumping around the converged solution for small errTol
    if iter > 2000 
        wtOldK      = 0.9900;
        wtOldC      = 0.9900;
        wtOldLambda = 0.9900;
    elseif iter > 5000
        wtOldK      = 0.9950;
        wtOldC      = 0.9950;
        wtOldLambda = 0.9950;
    end
end

if err2 <= errTol 
    fprintf("Model converged in %d iterations after %.2fs. MSE: %.12f\n", iter, timer, err2);
else
    fprintf("Model failed to converged after %d iterations in %.2fs. MSE: %.12f\n", iter, timer, err2);
end

%% Post Convergence (Solution) Analysis

saveas(gcf, './Figures/RTM_Soln_Path.png');

RegAnalysis = struct(); % store results 
vK_adj = vK(BurnT:end-BurnT); % only consider burn period 
for iA = 1:params.pNA     
    % create regression variables (lagged log K) 
    log_vK = log(vK_adj(find(ivA(BurnT:end-BurnT) == iA)));
    RegAnalysis.vK{iA} = log_vK(2:end);
    RegAnalysis.vK_Past{iA} = log_vK(1:end-1);
    % regression analysis conditioning on productivity 
    RegAnalysis.y{iA} = RegAnalysis.vK{iA}; 
    RegAnalysis.X{iA} = RegAnalysis.vK_Past{iA};
    RegAnalysis.X{iA} = [ones(size(RegAnalysis.X{iA})), RegAnalysis.X{iA}];
    [RegAnalysis.b{iA}, RegAnalysis.bint{iA}, RegAnalysis.r{iA}, RegAnalysis.rint{iA}, RegAnalysis.stats{iA}] = regress(RegAnalysis.y{iA}, RegAnalysis.X{iA});
    RegAnalysis.MSE{iA} = mean(RegAnalysis.r{iA}.^2);
    RegAnalysis.yhat{iA} = RegAnalysis.b{iA}(1) + RegAnalysis.b{iA}(2) .* RegAnalysis.X{iA}(:,2);
end

% figures 
subplot(3,3,1); 
scatter(RegAnalysis.vK_Past{1}, RegAnalysis.vK{1}, 6, 'blue','filled'); grid on;hold on; 
plot(RegAnalysis.X{1}(:,2), RegAnalysis.yhat{1}, 'r--');hold off; 
ylabel("Log(K')");xlabel('Log(K)');legend('A=1');
subplot(3,3,2);
scatter(RegAnalysis.vK_Past{2}, RegAnalysis.vK{2}, 6, 'blue','filled'); grid on;hold on; 
plot(RegAnalysis.X{2}(:,2), RegAnalysis.yhat{2}, 'r--');hold off; 
ylabel("Log(K')");xlabel('Log(K)');legend('A=2');
subplot(3,3,3);
scatter(RegAnalysis.vK_Past{3}, RegAnalysis.vK{3}, 6, 'blue','filled'); grid on;hold on; 
plot(RegAnalysis.X{3}(:,2), RegAnalysis.yhat{3}, 'r--');hold off; 
ylabel("Log(K')");xlabel('Log(K)');legend('A=3');
subplot(3,3,4);
scatter(RegAnalysis.vK_Past{4}, RegAnalysis.vK{4}, 6, 'blue','filled'); grid on;hold on; 
plot(RegAnalysis.X{4}(:,2), RegAnalysis.yhat{4}, 'r--');hold off; 
ylabel("Log(K')");xlabel('Log(K)');legend('A=4');
subplot(3,3,5);
scatter(RegAnalysis.vK_Past{5}, RegAnalysis.vK{5}, 6, 'blue','filled'); grid on;hold on; 
plot(RegAnalysis.X{5}(:,2), RegAnalysis.yhat{5}, 'r--');hold off; 
ylabel("Log(K')");xlabel('Log(K)');legend('A=5');
subplot(3,3,6);
scatter(RegAnalysis.vK_Past{6}, RegAnalysis.vK{6}, 6, 'blue','filled'); grid on;hold on; 
plot(RegAnalysis.X{6}(:,2), RegAnalysis.yhat{6}, 'r--');hold off; 
ylabel("Log(K')");xlabel('Log(K)');legend('A=6');
subplot(3,3,7);
scatter(RegAnalysis.vK_Past{7}, RegAnalysis.vK{7}, 6, 'blue','filled'); grid on;hold on; 
plot(RegAnalysis.X{7}(:,2), RegAnalysis.yhat{7}, 'r--');hold off; 
ylabel("Log(K')");xlabel('Log(K)');legend('A=7');
saveas(gcf, './Figures/K_scatter_plots.png'); 

% simulate K path -  actual nonlinear path vs log-linear prediction
size(vK)
vSimK = zeros(size(vK));
vSimK(1) = vK(1);
for t = 2:T 
    vSimK(t) = exp(RegAnalysis.b{ivA(t)}(1) + RegAnalysis.b{ivA(t)}(2)*log(vSimK(t-1)));
end
vDevK = abs(vK - vSimK);
maxDev = find(vDevK == max(vDevK));

figure;
plot(1:T, vK, 'b-', 'LineWidth',1); hold on; 
plot(1:T, vSimK, 'r-.', 'LineWidth',1); grid on;
yline(ss.K, 'k-');xlim([BurnT, T-BurnT]);ylabel('K');
xline(maxDev, 'k--');hold off;
legend('Actual', 'Log-lin Forecast', '','Max Deviation', 'Location', 'northeast');
saveas(gcf, './Figures/LogLinForecast.png');

% log variables 
vLogK = log(vK);
vLogI = log(vI);
vLogC = log(vC);
vLogL = log(vL);
vLogY = log(vY);

% relative variances 
std_Y2Y = std(vLogY)/std(vLogY);
std_K2Y = std(vLogK)/std(vLogY);
std_L2Y = std(vLogL)/std(vLogY);
std_I2Y = std(vLogI)/std(vLogY);
std_C2Y = std(vLogC)/std(vLogY);

% auto corr 
[acf, lags] = xcorr(vLogY - mean(vLogY), 'coeff');
acf = acf(lags>=0);
lags = lags(lags>=0);
autoCorr_Y = acf(2);

[acf, lags] = xcorr(vLogI - mean(vLogI), 'coeff');
acf = acf(lags>=0);
lags = lags(lags>=0);
autoCorr_I = acf(2);

[acf, lags] = xcorr(vLogC - mean(vLogC), 'coeff');
acf = acf(lags>=0);
lags = lags(lags>=0);
autoCorr_C = acf(2);

[acf, lags] = xcorr(vLogL - mean(vLogL), 'coeff');
acf = acf(lags>=0);
lags = lags(lags>=0);
autoCorr_L = acf(2);

[acf, lags] = xcorr(vLogK - mean(vLogK), 'coeff');
acf = acf(lags>=0);
lags = lags(lags>=0);
autoCorr_K = acf(2);

% co-movements 
corr_YC = corrcoef(vLogY, vLogC);
corr_YK = corrcoef(vLogY, vLogK);
corr_YI = corrcoef(vLogY, vLogI);
corr_YL = corrcoef(vLogY, vLogL);

% create table 
varNames = {'Y', 'I', 'C', 'L', 'K'};
rel_std = [std_Y2Y; std_I2Y; std_C2Y; std_L2Y; std_K2Y];
autoCorr = [autoCorr_Y; autoCorr_I; autoCorr_C; autoCorr_L; autoCorr_K];
corrCoef = [1.000; corr_YI(1,2); corr_YC(1,2); corr_YL(1,2); corr_YK(1,2)];
T = table(rel_std, autoCorr, corrCoef,...
    'RowNames', varNames,...
    'VariableNames', {'RelSd_to_Y','AutoCorr','CorrCoeff_to_Y'});
disp('Business cycle Statistics')
disp(T)