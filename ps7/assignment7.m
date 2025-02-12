%% Aiyagari PFI Approach 

%% Housekeeping 
clear;clc;close all;
addpath('./Figures')
addpath('./Functions')

%% Model Fundamentals 

% parameters 
p.pBeta = 0.990;
p.pAlpha = 0.330;
p.pDelta = 0.025;
p.pMu = 0.600;
p.pFrisch = 1.000;
p.pRiskAversion = 1.000;
p.pEta = 6.400;
p.pRhoA = 0.950;
p.pRhoZ = 0.900;
p.pSigmaA = 0.009;
p.pSigmaZ = 0.05;
p.pNA = 7;
p.pNz = 7;
p.paMin = 1.000;
p.paMax = 150.000;
p.pNa1 = 100;
p.pNa2 = 200;
p.pCurve = 5;

% steady state 
ss.A = 1;
ss.L = 1/3;

% unpack parameters - easy of code notation (struct is used in functions)
pBeta=p.pBeta;pAlpha=p.pAlpha;pDelta=p.pDelta;pMu=p.pMu;pFrisch=p.pFrisch;pRisk=p.pRiskAversion;pEta=p.pEta; 

% discretize TFP shock process and idiosyncratic shock process 
[vGridA, mPA] = fnTauchen(p.pRhoA,p.pSigmaA,0.0,p.pNA);
vGridA = exp(vGridA);
[vGridz, mPz] = fnTauchen(p.pRhoZ,p.pSigmaZ,0.0,p.pNz);
vGridz = exp(vGridz);

% set asset grid (asset-space) 
x1 = linspace(0,0.5,p.pNa1);
x2 = linspace(0,0.5,p.pNa2);
y1 = x1.^p.pCurve / max(x1.^p.pCurve);
y2 = x2.^p.pCurve / max(x2.^p.pCurve);
vGrida1 = p.paMin + (p.paMax - p.paMin).*y1; 
vGrida2 = p.paMin + (p.paMax - p.paMin).*y2;
ma = repmat(vGrida1', 1, p.pNz); % matrix of asset states over z-states
mz = repmat(vGridz', p.pNa1, 1); % matrix of productivity states over a-states %% CHECK THE DIMENSION HERE!!

%% solve for a SRCE using PFI Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fixed point method requires guessed solution and iterative convergence 
% guess (initial) solution 
K = 1.1776; % agg capital 
L = ss.L; % agg labour supply 
mLambda = zeros(p.pNa1, p.pNz); % state contingent savings constraint (lower bound)
mPolaprime = repmat(0.01 .* vGrida1', 1, p.pNz); % given pFrisch = pSigma = 1 we only need to guess a'(a,z) (initial) solution

% loop parameters 
wtOldK = 0.9900;
wtOldL = 0.9900;
wtOldmPolaprime = 0.9000;
wtOldLambda = 0.9000;
errTolSRCE = 1e-8;
errTolDist = 1e-8;
MaxIter = 10000;

% initialise 
Knew = 0;
mLambdanew = zeros(size(mLambda));
mPolaprimenew = zeros(size(mPolaprime));
mPolaprimenew2 = zeros(p.pNa2, p.pNz);
mCurrDist = ones(p.pNa1, p.pNz) ./ (p.pNa1*p.pNz);

% PFI Loop solve for SRCE
iter = 1;
err1 = 10;
tic;
while err1 > errTolSRCE && iter <= MaxIter
    
    % ==================================================
    % UPDATE EQUILIBRIUM PRICES
    % ==================================================

    r = pAlpha*ss.A*(K/L)^(pAlpha-1) - pDelta;
    w = (1-pAlpha).*ss.A.*(K/L)^(-pAlpha);
    
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
        mPolaprimeprime = interp1(vGrida1',mPolaprime(:,izprime),mPolaprime,'linear','extrap');
        
        % compute future auxillary variable over all current states (ia,iz)
        mMprime = ((1+rprime).*mPolaprime - (pMu/2).*((mPolaprimeprime-mPolaprime)./mPolaprime)./mPolaprime - mPolaprimeprime)./(wprime*zprime); %%%%%% CHECK EQ!!  

        % compute future optimal labour supply decision n' = n(ia'(ia,iz),iz') 
        mPolnprime = (-pEta.*mMprime + sqrt((pEta.*mMprime).^2 + 4*pEta))./(2*pEta); 

        % compute future optimal consumption decision c' = c(ia'(ia,iz),iz')
        mPolcprime = wprime*zprime./(pEta.*mPolnprime); %%%%%%%%%%%%%%%%%%%%%%%%% MAJOR - CHECK THIS EQ!!!!
        mPolcprime(mPolcprime<=1e-10) = 1e-10; %%%%%%%%%%%%%%%%%%%%%%%%% CONSTRAINT AT ZERO?? 

        % update beliefs 
        % compute expectation term cumulatively
        mMUinv = 1./mPolcprime;
        mBelief = mBelief +...
            pBeta .* repmat(mPz(:, izprime)', p.pNa1, 1).*...
            (1 + rprime - (pMu/2).*(mPolaprimeprime.^2 - mPolaprime.^2)./mPolaprime.^2).*mMUinv; %%%%%%%%%%%%%%% CHECK/SIMPLIFY THIS EQ AND SIGN!! 
    end
    
    % ==================================================
    % SOLVE FOR OPTIMAL POLICY RULES GIVEN BELIEFS 
    % ==================================================
    
    mPolc = (1 + pMu.*(mPolaprime-ma)./ma)./(mBelief + mLambda); %%%%%%%%%%%%% CHECK THIS EQUATION!!! 
    mPolc(mPolc<=1e-10) = 1e-10;
    mPoln = w.*mz./(pEta.*mPolc);
    
    % ==================================================
    % UPDATE N+1th ITERATION POLICY RULES AND CONSTRAINTS
    % ==================================================
    
    mLambdanew = (1+pMu.*(mPolaprime-ma)./ma)./((1+r).*ma + w.*mz.*mPoln - ma.*(pMu/2).*((mPolaprime-ma)./ma).^2 - mPolaprime) - mBelief; %%% CHECK THIS UPDATING!! 
    mPolaprimenew = (1+r).*ma + w.*mz.*mPoln - ma.*(pMu/2).*((mPolaprime-ma)./ma).^2 - mPolc;

    % impose friction (lower bound on savings constraint)
    mPolaprimenew(mPolaprimenew<=1) = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RELATE CONSTRAINTS TO 
    mLambdanew(mPolaprimenew>1) = 0;
    
    % ==================================================
    % GIVEN POLICY RULES COMPUTE STATIONARY DISTRIBUTION
    % ==================================================

    % interpolate policy rule on finer grid 
    % for iz = 1:p.pNz
    %     mPolaprimenew2(:, iz) = interp1(vGrida1, mPolaprimenew(:, iz), vGrida000000000000000000000000000000000000000000000000002, "linear", "extrap");
    % end

    % non-stochastic iterative histogram method
    iter2 = 1;
    err2 = 10;
    while err2 > errTolDist

        mNewDist = zeros(size(mCurrDist));

        for iz = 1:p.pNz 
            for ia = 1:p.pNa1
                
                % for state (ia,iz) interpolate optimal savings on asset grid
                aprime = mPolaprimenew(ia,iz);
                nLow = sum(vGrida1 <= aprime);
                nLow(nLow<=1) = 1;
                nLow(nLow>=length(vGrida1)) = length(vGrida1) - 1;
                nHigh = nLow + 1;
                wtLow = (vGrida1(nHigh)-aprime)/(vGrida1(nHigh)-vGrida1(nLow));
                wtLow(wtLow>1) = 1;
                wtLow(wtLow<0) = 0;
                
                % get current mass over state (ia,iz) 
                mass = mCurrDist(ia,iz);
                
                % update mass cumutatively for all possible future iz states
                for izprime = 1:p.pNz

                    mNewDist(nLow, izprime) = mNewDist(nLow, izprime) + ...
                        1 * mass * mPz(iz,izprime) * wtLow;
                    mNewDist(nHigh, izprime) = mNewDist(nHigh, izprime) + ...
                        1 * mass * mPz(iz,izprime) * (1-wtLow);
                end
            end
        end
        
        % compute error and update distribution 
        err2 = max(max(abs(mNewDist-mCurrDist)));
        mCurrDist = mNewDist;
        % if mod(iter2, 1)==0
        %     fprintf('Iter %d. Distribution Error: %.10f\n', iter2, err2)
        % end
        % iter2=iter2+1;
    end

    % ==================================================
    % AGGREGATION AND CLEAR MARKETS 
    % ==================================================
    
    vMargDista = sum(mCurrDist,2);
    Knew = vGrida1 * vMargDista;
    Lnew = sum(repmat(vGridz', p.pNa1, 1) .* mPoln .* mCurrDist, 'all');
    
    % ==================================================
    % UPDATING 
    % ==================================================
    
    % compute error
    errK = abs(K-Knew);
    errL = abs(L-Lnew);
    meanLambda = sum(mLambda .* mCurrDist, 'all');
    meanLambdaNew = sum(mLambdanew .* mCurrDist, 'all');
    errPolaprime = max(max(abs(mPolaprime - mPolaprimenew)));
    errLambda = max(max(abs(mLambda - mLambdanew)));
    
    err1 = mean([...
        errK; ...
        errL; ...
        abs(meanLambda - meanLambdaNew)]);

    % convex updating (smooth)
    K = wtOldK*K + (1-wtOldK)*Knew;
    L = wtOldL*L + (1-wtOldL)*Lnew;
    mPolaprime = wtOldmPolaprime.*mPolaprime + (1-wtOldmPolaprime).*mPolaprimenew;
    mLambda = wtOldLambda.*mLambda + (1-wtOldLambda).*mLambdanew;
    
    timer = toc;
    if mod(iter, 2) == 0
        
        fprintf('Iteration %d in %.2fs. Error: %.10f\n', iter, timer, err1);
        fprintf('errK: %.6f. errL: %.6f. errPolaprime: %.6f. errLambda: %.6f\n', errK, errL, errPolaprime, errLambda);
    end
    iter = iter+1;
end