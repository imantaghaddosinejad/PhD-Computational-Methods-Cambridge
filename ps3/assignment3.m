%% Assignment 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping
clear;
clc;

cd 'C:\Users\it315\Documents\phd_computational_macro\ps3';
addpath('./../ReplicationFiles/Aiyagari1994/Matlab/')
addpath('./../ReplicationFiles/Aiyagari1994/Matlab/Functions/')
addpath('./../ReplicationFiles/Aiyagari1994/Matlab/Figures/')
addpath('./Figures')

%% solve Aiyagari 1994
Main_VFI;

%% Ex. (f)



%% Ex. (g)

% creating the joint-transition matrix over states (a,z)
mJointTrans = zeros(pNa2*pNz, pNa2*pNz);
for ia = 1:pNa2
    for iz = 1:pNz
        iCurrState = (iz - 1)*pNa2 + ia; % corresponding i-th row for state (ia, iz)
        aprime = mPolaprime2(ia, iz); % optimal savings given current state (ia, iz)
        [LB, UB, wtLB, wtUB] = fnInterp1dGrid(aprime, vGrida2, pNa2); % interpolate aprime on asset grid
        
        for izfuture = 1:pNz
            iFutureStateLB = (izfuture - 1)*pNa2 + LB; % index for LB - future state 
            iFutureStateUB = (izfuture - 1)*pNa2 + UB; % index for UB - future state 
            mJointTrans(iCurrState, iFutureStateLB) = 1 * mPz(iz, izfuture) * wtLB; % prob. transition to LB
            mJointTrans(iCurrState, iFutureStateUB) = 1 * mPz(iz, izfuture) * wtUB; % prob. transition to UB
        end
    end
end

% simulate 10,000 households for 200 periods 
T = 200;
Nh = 10000;

mPathAssets = zeros(Nh, T+1); % placeholder for simulated path of assets 
mPathZ = zeros(Nh, T+1); % placeholder for simulated path of productivity
miPath = zeros(Nh, T+1); % placeholder for simulated paths joint-index over (a,z) 
miPath(:, 1) = randsample(pNa2*pNz, Nh, true, mCurrDist(:)); % set initial states for households 

for t = 1:T
    for iH = 1:Nh
        iHCurrState = miPath(iH, t); % houehold current state 
        iHFutureState = randsample(pNa2*pNz, 1, true, mJointTrans(iHCurrState, :)); % draw future state from joint-transition matrix 
        miPath(iH, t+1) = iHFutureState; % update path matrix for current household 
    end
end

% unpack (convert back) state indices to state values for each household over time
for t= 1:T+1
    [iA, iZ] = ind2sub([pNa2, pNz], miPath(:, t)); % convert index to (iA, iZ) state 
    mPathAssets(:, t) = vGrida2(iA); % update asset path
    mPathZ(:,t) = vGridz(iZ); % update productivity path
end

% plotting aggregate wealth path over time 
vPathAggAssets = sum(mPathAssets, 1);
figure;
plot(1:T, vPathAggAssets(1:T), "LineWidth", 1.5)
legend('Aggregate Wealth')
grid on;
saveas(gcf, "./Figures/agg_wealth.png")

%% Ex. (h) 

% find index of households in the bottom (top) 10% of wealth in T=100
Tpoint = 100;
Tstart = Tpoint + 1;
bottom_10pct_indx = find(mPathAssets(:, Tpoint) <= prctile(mPathAssets(:, Tpoint), 10));
top_10pct_indx = find(mPathAssets(:, Tpoint) >= prctile(mPathAssets(:, Tpoint), 90));

% find mean wealth (per year) of bottom (top) 10% households 
mean_bottom_10pct = mean(mPathAssets(bottom_10pct_indx, Tstart:T), 1)';
mean_top_10pct = mean(mPathAssets(top_10pct_indx, Tstart:T), 1)';

% compute standard error for bottom (top) 10% wealth (across all years) of households 
std_bottom_10pct = std(mPathAssets(bottom_10pct_indx, Tstart:T), 0, 2);
std_bottom_mean = mean(std_bottom_10pct);
std_top_10pct = std(mPathAssets(top_10pct_indx, Tstart:T), 0, 2);
std_top_mean = mean(std_top_10pct);

% confidence interval computation
ci_bottom_10pct = [mean_bottom_10pct - 1.96 .* std_bottom_mean, mean_bottom_10pct + 1.96 .* std_bottom_mean];
ci_top_10pct = [mean_top_10pct - 1.96 .* std_top_mean, mean_top_10pct + 1.96 .* std_top_mean];

% prepare variables
xtime = Tstart:T;
xtime = xtime';
bottom_ci_LB = ci_bottom_10pct(:, 1);
bottom_ci_UB = ci_bottom_10pct(:, 2);
top_ci_LB = ci_top_10pct(:, 1);
top_ci_UB = ci_top_10pct(:, 2);

% plot mean wealth and standard errors 
figure;
hold on;
fill([xtime; flipud(xtime)], [top_ci_LB; flipud(top_ci_UB)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(xtime, mean_top_10pct, 'b', 'LineWidth', 1.5);  
fill([xtime; flipud(xtime)], [bottom_ci_LB; flipud(bottom_ci_UB)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(xtime, mean_bottom_10pct, 'r', 'LineWidth', 1.5);  
hold off;
grid on;
legend('95% CI', 'Top 10%', '95% CI', 'Bottom 10%')
title('Mean Wealth')
saveas(gcf, "./Figures/mean_wealth.png")

%% Ex. (j)

% plotting individual wealth paths over time
meanAssets = mean(mean(mPathAssets)); 
[~, idx1] = min(abs(mPathAssets(:, 1) - 49));
[~, idx2] = min(abs(mPathAssets(:, 1) - 7));
h1 = [];
h2 = [];
dt_time = 51:100;

figure;
hold on;
for i = 1:Nh
    if i == idx1
        h1 = plot(dt_time, mPathAssets(idx1, dt_time), "LineWidth", 4, "color", [0.6, 0, 0]);
    elseif i == idx2
        continue;
    else
        plot(dt_time, mPathAssets(i, dt_time), "LineWidth", 0.35)        
    end
end
m1 = yline(meanAssets, "LineWidth", 4, "Color", [0, 0, 0]);
h2 = plot(dt_time, mPathAssets(idx2, dt_time), "LineWidth", 4, "Color", [0, 0, 0.7]);

legend([h1, h2, m1], {sprintf('Household %d', idx1), ...
    sprintf('Household %d',  idx2), ...
    sprintf('Mean %.2f', meanAssets)});
hold off;
grid on;
saveas(gcf, "./Figures/individual_wealth.png")

%% Ex. (j) 

% simulate 10,000 households for 100 periods 
T = 100;
Nh = 10000;

mPathAssets = zeros(Nh, T+1); % placeholder for simulated path of assets 
mPathZ = zeros(Nh, T+1); % placeholder for simulated path of productivity
miPath = zeros(Nh, T+1); % placeholder for simulated paths joint-index over (a,z) 
miPath(:, 1) = randsample(pNa2*pNz, Nh, true, mCurrDist(:)); % set initial states for households 

for t = 1:T
    for iH = 1:Nh
        iHCurrState = miPath(iH, t); % houehold current state 
        iHFutureState = randsample(pNa2*pNz, 1, true, mJointTrans(iHCurrState, :)); % draw future state from joint-transition matrix 
        miPath(iH, t+1) = iHFutureState; % update path matrix for current household 
    end
end

% unpack (convert back) state indices to state values for each household over time
for t= 1:T+1
    [iA, iZ] = ind2sub([pNa2, pNz], miPath(:, t)); % convert index to (iA, iZ) state 
    mPathAssets(:, t) = vGrida2(iA); % update asset path
    mPathZ(:,t) = vGridz(iZ); % update productivity path
end

poolwealth = mPathAssets(:, 51:100);
poolwealth = poolwealth(:);
poolwealth = sort(poolwealth);
cum_poolwealth = cumsum(poolwealth);
total_poolwealth = sum(poolwealth);
wealth_share = cum_poolwealth ./ total_poolwealth;

n = length(poolwealth);
x_lorenz = (1:n) ./ n ; 

x_lorenz = [0, x_lorenz];
y_lorenz = [0, wealth_share'];

gini = 1 - 2 * trapz(linspace(0, 1, length(y_lorenz)), y_lorenz);

plot(x_lorenz, y_lorenz, 'b-', 'LineWidth', 2);
hold on;
plot([0, 1], [0, 1], 'r--', 'LineWidth', 1);
hold off;
grid on;
xlabel('Cumulative Populaiton Proportion');
ylabel('Cumulative Wealth Proportion');
legend('Lorenz Curve (TFP=1', 'Line of Equality')
axis([0 1 0 1]);

text(0.8, 0.1, ['Gini Coefficient: ', num2str(gini, '%.3f')], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', [1 1 0.9], ...  % Light yellow background
    'EdgeColor', [0.6 0.6 0.6], ...    % Gray border
    'FontSize', 10);
saveas(gcf, "./Figures/lorenz_curve.png")

%% Ex. 

mJointTrans = zeros(pNa2*pNz, pNa2*pNz);
for ia = 1:pNa2
    for iz = 1:pNz
        iCurrState = (iz - 1)*pNa2 + ia; % corresponding i-th row for state (ia, iz)
        aprime = mPolaprime2(ia, iz); % optimal savings given current state (ia, iz)
        [LB, UB, wtLB, wtUB] = fnInterp1dGrid(aprime, vGrida2, pNa2); % interpolate aprime on asset grid
        
        for izfuture = 1:pNz
            iFutureStateLB = (izfuture - 1)*pNa2 + LB; % index for LB - future state 
            iFutureStateUB = (izfuture - 1)*pNa2 + UB; % index for UB - future state 
            mJointTrans(iCurrState, iFutureStateLB) = 1 * mPz(iz, izfuture) * wtLB; % prob. transition to LB
            mJointTrans(iCurrState, iFutureStateUB) = 1 * mPz(iz, izfuture) * wtUB; % prob. transition to UB
        end
    end
end

% simulate 10,000 households for 100 periods 
T = 100;
Nh = 10000;

mPathAssets = zeros(Nh, T+1); % placeholder for simulated path of assets 
mPathZ = zeros(Nh, T+1); % placeholder for simulated path of productivity
miPath = zeros(Nh, T+1); % placeholder for simulated paths joint-index over (a,z) 
miPath(:, 1) = randsample(pNa2*pNz, Nh, true, mCurrDist(:)); % set initial states for households 

for t = 1:T
    for iH = 1:Nh
        iHCurrState = miPath(iH, t); % houehold current state 
        iHFutureState = randsample(pNa2*pNz, 1, true, mJointTrans(iHCurrState, :)); % draw future state from joint-transition matrix 
        miPath(iH, t+1) = iHFutureState; % update path matrix for current household 
    end
end

% unpack (convert back) state indices to state values for each household over time
for t= 1:T+1
    [iA, iZ] = ind2sub([pNa2, pNz], miPath(:, t)); % convert index to (iA, iZ) state 
    mPathAssets(:, t) = vGrida2(iA); % update asset path
    mPathZ(:,t) = vGridz(iZ); % update productivity path
end

poolwealth = mPathAssets(:, 51:100);

poolwealth = poolwealth(:);
poolwealth = sort(poolwealth);
cum_poolwealth = cumsum(poolwealth);
total_poolwealth = sum(poolwealth);
wealth_share = cum_poolwealth ./ total_poolwealth;

n = length(poolwealth);
x_lorenz = (1:n) ./ n ; 

x_lorenz = [0, x_lorenz];
y_lorenz = [0, wealth_share'];

gini = 1 - 2 * trapz(linspace(0, 1, length(y_lorenz)), y_lorenz);

plot(x_lorenz, y_lorenz, 'b-', 'LineWidth', 2);
hold on;
plot([0, 1], [0, 1], 'r--', 'LineWidth', 1);
hold off;
grid on;
xlabel('Cumulative Populaiton Proportion');
ylabel('Cumulative Wealth Proportion');
legend('Lorenz Curve (TFP = 1.2)', 'Line of Equality')
axis([0 1 0 1]);

text(0.8, 0.1, ['Gini Coefficient: ', num2str(gini, '%.3f')], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', [1 1 0.9], ...  % Light yellow background
    'EdgeColor', [0.6 0.6 0.6], ...    % Gray border
    'FontSize', 10);
saveas(gcf, "./Figures/lorenz_curve2.png")