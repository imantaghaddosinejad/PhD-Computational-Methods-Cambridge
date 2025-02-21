%% Aiyagari PFI Approach 

%% Housekeeping 
close all;clc;
clear variables;
addpath('./Figures')
addpath('./Functions')

%% Model Fundamentals 

% parameters 
p.pBeta = 0.96;
p.pAlpha = 0.33;
p.pDelta = 0.08;
p.pMu = 0;%0.600;
p.pFrisch = 1.0;
p.pRiskAversion = 1.0;
p.pEta = 7.60;
p.pRhoA = 0.950;
p.pRhoZ = 0.90;
p.pSigmaA = 0.009;
p.pSigmaZ = 0.05;
p.pNA = 7;
p.pNz = 7;
p.paMin = 1.0;
p.paMax = 200;
p.pNa = 100;
p.pCurve = 7;

% steady state 
ss.A = 1;
ss.L = .33;

% unpack parameters - easy of code notation (struct is used in functions)
pBeta=p.pBeta;pAlpha=p.pAlpha;pDelta=p.pDelta;pMu=p.pMu;pFrisch=p.pFrisch;pRisk=p.pRiskAversion;pEta=p.pEta; 

% discretize TFP shock process 
[vGridA, mPA] = fnTauchen(p.pRhoA,p.pSigmaA,0.0,p.pNA);
vGridA = exp(vGridA);

% discretize idiosyncratic shock process 
[vGridz, mPz] = fnTauchen(p.pRhoZ,p.pSigmaZ,0.0,p.pNz);
vGridz = exp(vGridz);

% set asset grid (asset-space) 
x = linspace(0,0.5,p.pNa);
y = x.^p.pCurve / max(x.^p.pCurve);
vGrida = p.paMin + (p.paMax - p.paMin).*y; 

% auxillary matrices 
ma = repmat(vGrida', 1, p.pNz); % matrix of asset states over z-states
mz = repmat(vGridz', p.pNa, 1); % matrix of productivity states over a-states 

%% solve for a SRCE using PFI Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fixed point method requires guessed solution and iterative convergence 
% guess (initial) solution
K = 3.5; % agg capital 
L = ss.L; % agg labour supply 
%mPolaprime = repmat(0.01.* vGrida', 1, p.pNz); % savings policy rule 
mPolaprime = ones(p.pNa, p.pNz);
mPoln = zeros(size(mPolaprime)); % labour supply policy rule 

% equilibrium objects 
%rinit = pAlpha*ss.A*(K/L)^(pAlpha-1)-pDelta;
%winit = (1-pAlpha)*ss.A*(K/L)^pAlpha;
%mPolc = (1+rinit).*ma + winit.*mz.*mPoln - mPolaprime - (pMu/2).*ma.*((mPolaprime-ma)./ma).^2;

mLambda = zeros(p.pNa, p.pNz); % state contingent lagrangian multiplier for NNC 
%mPolc = repmat(0.01.* vGrida', 1, p.pNz); % policy rule for consumption
mPolc = ones(p.pNa,p.pNz);

% loop parameters 
wtOldK = 0.9000;
wtOldL = 0.9000;
wtOldmPolaprime = 0.9000;
wtOldLambda = 0.9000;
wtOldmPolc = 0.9000;
errTolSRCE = 1e-8;
errTolDist = 1e-10;
MaxIter = 50000;

% initialise 
Knew = 0;
mLambdanew = zeros(size(mLambda));
mPolaprimenew = zeros(size(mPolaprime));
mCurrDist = ones(p.pNa,p.pNz)./(p.pNa*p.pNz);

% PFI Loop solve for SRCE
iter = 1;
err1 = 10;
printinterval=35;
tic;
while err1 > errTolSRCE && iter <= MaxIter
    
    % ==================================================
    % UPDATE EQUILIBRIUM PRICES
    % ==================================================

    r = pAlpha*ss.A*(K/L)^(pAlpha-1) - pDelta;
    w = (1-pAlpha).*ss.A.*(K/L)^(pAlpha);
    
    % ==================================================
    % SOLVE BACKWARDS FOR BELIEFS GIVEN Nth POLICY RULES
    % ==================================================

    % compute belief over states (ia,iz) 
    mBelief = 0; % reset expectation term 
    for izprime = 1:p.pNz
        
        % future prices 
        rprime = r;
        wprime = w;

        % future (realised) state 
        zprime = vGridz(izprime);

        % interpolate policy function to compute a'' = a''(ia',iz') over all current states (ia,iz)
        mPolaprimeprime = interp1(vGrida', squeeze(mPolaprime(:,izprime)), squeeze(mPolaprime),"linear","extrap");
        
        % cost term 
        mPsi1 = (pMu/2).*((mPolaprimeprime./mPolaprime).^2-1);

        % compute future auxillary variable over all current states (ia,iz)
        mMprime = ((1+rprime).*mPolaprime - mPolaprimeprime - (pMu/2).*((mPolaprimeprime-mPolaprime)./mPolaprime).^2.*mPolaprime)./(wprime*zprime); 

        % compute future optimal labour supply decision n' = n(ia'(ia,iz),iz') 
        mPolnprime = (-pEta.*mMprime + sqrt((pEta.*mMprime).^2 + 4*pEta))./(2*pEta); 

        % compute future optimal consumption decision c' = c(ia'(ia,iz),iz')
        %mPolcprime = (wprime.*zprime)./(pEta.*mPolnprime); 
        mPolcprime =  wprime.*zprime*mPolnprime + (1+rprime).*squeeze(mPolaprime) - mPolaprimeprime  ...
           - (pMu/2).*((mPolaprimeprime-mPolaprime)./mPolaprime).^2.*mPolaprime;
        mPolcprime(mPolcprime<=0) = 1e-10;

        % update beliefs 
        % compute expectation term cumulatively
        mMUinv = 1./mPolcprime;
        mBelief = mBelief + repmat(mPz(:, izprime)',p.pNa,1).*(1+rprime+mPsi1).*mMUinv; 
    end
    
    % ==================================================
    % SOLVE FOR OPTIMAL POLICY RULES GIVEN BELIEFS 
    % ==================================================
    
    mBelief = pBeta.*mBelief;
    c = (1+pMu.*(mPolaprime-ma)./ma)./(mBelief+mLambda);
    c(c<=0) = 1e-10;
    mPoln = (w.*mz)./(pEta.*c);
    
    % ==================================================
    % UPDATE N+1th ITERATION POLICY RULES AND CONSTRAINTS
    % ==================================================
    
    mLambdanew = (1+pMu.*(mPolaprime-ma)./ma)./(mPolc) - mBelief;
    %mLambdanew = (1+pMu.*(mPolaprime-ma)./ma)./(w.*mz.*mPoln + (1+r).*ma - mPolaprime - (pMu/2).*((mPolaprime-ma)./ma).^2.*ma) - mBelief;
    mPolaprimenew =  w.*mz.*mPoln + (1+r).*ma - (pMu/2).*((mPolaprime-ma)./ma).^2.*ma - c;

    % impose friction (lower bound on savings constraint)
    mPolaprimenew(mPolaprimenew<=p.paMin) = p.paMin;
    mLambdanew(mPolaprimenew>p.paMin) = 0;
    %mPolc = c;
    
    % ==================================================
    % GIVEN POLICY RULES COMPUTE STATIONARY DISTRIBUTION
    % ==================================================

    % non-stochastic iterative histogram method
    iter2 = 1;
    err2 = 10;
    while err2 > errTolDist

        mNewDist = zeros(size(mCurrDist));
        for iz = 1:p.pNz 
            for ia = 1:p.pNa
                
                % for state (ia,iz) interpolate optimal savings on asset grid
                aprime = mPolaprimenew(ia,iz);
                nLow = sum(vGrida < aprime);
                nLow(nLow<=1) = 1;
                nLow(nLow>=length(vGrida)) = length(vGrida) - 1;
                nHigh = nLow + 1;
                wtLow = (vGrida(nHigh)-aprime)/(vGrida(nHigh)-vGrida(nLow));
                wtLow(wtLow>1) = 1;
                wtLow(wtLow<0) = 0;
                wtHigh = 1-wtLow;
                
                % get current mass over state (ia,iz) 
                mass = mCurrDist(ia,iz);
                
                % update mass cumutatively for all possible future iz states
                for izprime = 1:p.pNz

                    mNewDist(nLow, izprime) = mNewDist(nLow, izprime) + ...
                        1 * mass * mPz(iz,izprime) * wtLow;
                    
                    mNewDist(nHigh, izprime) = mNewDist(nHigh, izprime) + ...
                        1 * mass * mPz(iz,izprime) * (wtHigh);
                end
            end
        end
        
        % compute error and update distribution 
        err2 = max(abs(mNewDist-mCurrDist),[],'all');
        mCurrDist = mNewDist;

        % if mod(iter2, 25) == 0
        %      fprintf('dist error %.10f\n', err2);
        % end

        iter2 = iter2+1;
    end

    % fprintf('Dist Converged: %d iterations\n', iter2)

    % ==================================================
    % AGGREGATION AND CLEAR MARKETS 
    % ==================================================
    
    % market clearing
    vMargDista          = sum(mCurrDist,2);
    Knew                = vGrida*vMargDista;
    Lnew                = sum(mz.*mPoln.*mCurrDist,'all');
    avg_mLambda         = sum(mLambda.*mCurrDist,'all');
    avg_mLambdanew      = sum(mLambdanew.*mCurrDist,'all');
    avg_mPolaprime      = sum(mPolaprime.*mCurrDist,'all');
    avg_mPolaprimenew   = sum(mPolaprimenew.*mCurrDist,'all');
    avg_mPolc           = sum(mPolc.*mCurrDist, 'all');
    avg_mPolcnew        = sum(c.*mCurrDist, 'all');

    % compute error 
    err1 = mean(abs([...
        Knew                - K;...
        Lnew                - L;...
        avg_mLambdanew      - avg_mLambda;...
        avg_mPolaprimenew   - avg_mPolaprime;...
        avg_mPolcnew        - avg_mPolc]), 'all');

    % ==================================================
    % UPDATING 
    % ==================================================

    % convex updating (smooth)
    K           = wtOldK*K                      + (1-wtOldK)*Knew;
    L           = wtOldL*L                      + (1-wtOldL)*Lnew;
    mPolaprime  = wtOldmPolaprime.*mPolaprime   + (1-wtOldmPolaprime).*mPolaprimenew;
    mLambda     = wtOldLambda.*mLambda          + (1-wtOldLambda).*mLambdanew;
    mPolc       = wtOldmPolc.*mPolc             + (1-wtOldmPolc).*c;

    % ==================================================
    % PRINT PROGRESS
    % ==================================================

    timer = toc;
    if mod(iter, printinterval) == 0
        fprintf('-------------------------------------------------------\n');
        fprintf('market clearing results: iteration %d in %.2fs\n', iter, timer);
        fprintf('max error: %.12f \n', err1);
        fprintf('r: %.6f    w: %.6f     K: %.6f    L: %.6f\n', r, w, K, L);
        fprintf('err_lambda: %.8f   min_lambda: %.4f    Dist Converged: %d iters\n', abs(avg_mLambda-avg_mLambdanew), min(min(mLambda)), iter2)
    end

    % if mod(iter, 500) == 0
    % plot(vGrida, sum(mCurrDist,2),'LineWidth',1); grid on;drawnow;
    % end
    iter = iter+1;
end

%%
%====================
% SOLUTION ANALYSIS 
%====================
% compute VF 
mU = log(mPolc) - pEta.*mPoln.^(1+1/pFrisch)./(1+1/pFrisch);
mV = zeros(size(mU));
err = 10; 
iter = 1;
while err > 1e-6 
    mExpV = mV * mPz';
    mVnew = mU + pBeta.*mExpV;
    err = max(max(abs(mVnew-mV)));
    mV = mVnew;
end
fprintf('VF converged after %d iterations!\n', iter)

% save file 

%%
% plots

% 1 
subplot(1,2,1);
plot(vGrida,mCurrDist,'LineWidth',1.5); grid on;xlabel('a');ylabel('f(a)');
legend('z_1','z_2','z_3','z_4','z_5','z_6','z_7',"FontSize",20);
subplot(1,2,2);
plot(vGrida,mPolaprime(:,1)); hold on; grid on;
plot(vGrida,mPolaprime(:,end));ylabel("a'"); xlabel("a")
legend('Lowest z','Highest z','location','southeast',"FontSize",20);
hold off;
%saveas(gcf,'../figures/spec2/policy.jpg');

% 2
figure;
plot(vGrida, sum(mCurrDist,2),'LineWidth',1.5);grid on;
%saveas(gcf,'../figures/spec2/marg_dist_a.jpg');

% 3
figure;
subplot(1,2,1);
plot(vGrida, mPolc(:,1)); hold on; grid on;
plot(vGrida, mPolc(:,3));
plot(vGrida, mPolc(:,7));
legend('z1','z3','z7','Location','southeast',"FontSize",15);ylabel('c');xlabel('a');
hold off; 
subplot(1,2,2);
plot(vGrida, mPoln(:,1)); hold on; grid on;
plot(vGrida, mPoln(:,3));
plot(vGrida, mPoln(:,7));
legend('z1','z3','z7','Location','northeast',"FontSize",15);ylabel('n');xlabel('a');
%saveas(gcf,'../figures/spec2/cons_laborsupply_pol.jpg');

% 4
figure;
subplot(2,1,1);
plot(vGrida, mU(:,1)); grid on; hold on;
plot(vGrida, mU(:,2));
plot(vGrida, mU(:,5));
plot(vGrida, mU(:,6));
plot(vGrida, mU(:,7)); hold off;
legend('z1','z2','z5','z6','z7','Location','southeast',"FontSize",10);ylabel('U');xlabel('a');
subplot(2,1,2)
plot(vGrida, mV(:,1)); grid on; xlabel('a');ylabel('V'); hold on;
plot(vGrida, mV(:,2));
plot(vGrida, mV(:,5));
plot(vGrida, mV(:,6));
plot(vGrida, mV(:,7)); hold off;
legend('z1','z2','z5','z6','z7','Location','southeast',"FontSize",10);
%saveas(gcf,'../figures/spec2/utility_VF.jpg');




%%
figure;
subplot(1,2,1);
plot(vGrida, soln_spec2_hl.mpolc(:,1)); hold on; grid on;
plot(vGrida, soln_spec2_hl.mpolc(:,3));
plot(vGrida, soln_spec2_hl.mpolc(:,7));
legend('z1','z3','z7','Location','southeast',"FontSize",15);ylabel('c');xlabel('a');
hold off; 
subplot(1,2,2);
plot(vGrida, soln_spec2_hl.mpoln(:,1)); hold on; grid on;
plot(vGrida, soln_spec2_hl.mpoln(:,3));
plot(vGrida, soln_spec2_hl.mpoln(:,7));
legend('z1','z3','z7','Location','northeast',"FontSize",15);ylabel('n');xlabel('a');
%saveas(gcf,'../figures/spec2/cons_laborsupply_pol.jpg');


soln_hl = load('./Hanbeak_Lee/soln_ss_spec2.mat');

fprintf('-------------------------------------------------------\n');
fprintf('market clearing results: iteration %d in %.2fs\n', iter, timer);
fprintf('max error: %.12f \n', err1);
fprintf('r: %.6f    w: %.6f     K: %.6f    L: %.6f\n', soln_hl.r, soln_hl.w, soln_hl.K, soln_hl.supplyL);
fprintf('err_lambda: %.8f   min_lambda: %.4f\n', abs(avg_mLambda-avg_mLambdanew), min(min(mLambda)))


