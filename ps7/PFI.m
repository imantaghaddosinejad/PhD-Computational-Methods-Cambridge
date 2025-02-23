%% Housekeeping 
close all;clc;
clear variables;
% addpath('./Figures')
addpath('./Functions')

%% Model Fundamentals 
% model parameters 
p.Beta         = 0.99;
p.Alpha        = 0.36;
p.Delta        = 0.025;
p.Theta        = 0.5;
p.Frisch       = 1.0;
p.Eta          = 7.6;
p.RiskAversion = 1.0;
p.Mu           = 0.0;

% other parameters
p.RhoA     = 0.95;
p.Rhoz      = 0.90;
p.SigmaA   = 0.009;
p.Sigmaz   = 0.05;
p.NA       = 7;
p.Nz       = 7;
p.Mina      = 1e-10;
p.Maxa      = 200;
p.Na       = 100;
p.Curve    = 7;

% unpack params 
pbeta=p.Beta;palpha=p.Alpha;pdelta=p.Delta;ptheta=p.Theta;pfrisch=p.Frisch;prisk=p.RiskAversion;peta=p.Eta;pmu=p.Mu;

%% Grids and Random Processes 
% asset grid 
x = linspace(0,0.5,p.Na);
y = x.^p.Curve / max(x.^p.Curve);
vGrida = p.Mina + (p.Maxa - p.Mina).*y;

% idiosyncratic productivity grid
[vGridz, mPz] = fnTauchen(p.Rhoz,p.Sigmaz,0.0,p.Nz);
vGridz = exp(vGridz);

% aggregate productivity grid
[vGridA, mPA] = fnTauchen(p.RhoA,p.SigmaA,0.0,p.NA);
vGridA = exp(vGridA);

%% SRCE
% equilibrium objects 
mPolc           = repmat(0.01.*vGrida', 1, p.Nz);
K_new           = 0;
L_new           = 0;
mPolaprime_new  = zeros(size(mPolc)); 
mPoln_new       = zeros(size(mPolc));
mLambda_new     = zeros(size(mPolc));
A               = 1.0;

%  initial guess
K           = 11;
L           = 0.28;
mPolaprime  = ones(size(mPolc)); 
mPoln       = ones(size(mPolc)).*L;
mLambda     = zeros(size(mPolc));
mCurrDist   = ones(p.Na,p.Nz)/(p.Na*p.Nz);

% auxillary objects 
mgrida = repmat(vGrida', 1, p.Nz); 
mgridz = repmat(vGridz', p.Na, 1);

% loop parameters 
tolGE = 1e-8;
tolDist = 1e-10;
wtOld1 = 0.8000;
wtOld2 = 0.8000;
wtOld3 = 0.8000;
wtOld4 = 0.8000;
wtOld5 = 0.8000;

%=========================
% GE LOOP
%=========================

iter1 = 1;
err1 = 10;
dispiter = 25;
tic;
while err1 > tolGE 

    %=========================
    % UPDATE PRICES
    %=========================

    r = palpha * A * (K/L)^(palpha-1) - pdelta;
    w = (1-palpha) * A * (K/L)^palpha;
    mmu = r+pdelta;

    %=========================
    % UPDATE BELIEFS 
    %=========================

    mExp = 0;
    for izprime = 1:p.Nz 
        
        % future realised state 
        rprime = r;
        wprime = w;
        zprime = vGridz(izprime);
        
        % interpolate policy rules  
        aprimeprime = interp1(vGrida', squeeze(mPolaprime(:,izprime)), mPolaprime, 'linear', 'extrap');
        nprime = interp1(vGrida', squeeze(mPoln(:,izprime)), mPolaprime, 'linear', 'extrap');
        cprime = wprime.*zprime.*nprime + (1+rprime).*mPolaprime - aprimeprime;
        cprime(cprime<=0) = 1e-10;
        muprime = 1./cprime.^prisk;
        mExp = mExp + repmat(mPz(:,izprime)', p.Na, 1).*(1+rprime).*muprime;
    end

    %==========================
    % OPTIMAL POLICY GIVEN BELIEFS
    %==========================

    mExp = pbeta*mExp;
    c = (1./(mExp+mLambda)).^(1/prisk);
    c(c<=0) = 1e-10;
    n = ((w*mgridz)./(peta*c.^prisk)).^pfrisch;
    n(n>=1) = 1;
    n(n<=0) = 0;

    %==========================
    % UPDATE POLICY RULES
    %==========================

    mLambda_new = 1./((1+r)*mgrida + w*mgridz.*n - mPolaprime).^prisk - mExp;
    mPolaprime_new = (1+r)*mgrida + w*mgridz.*n - c;
    mPoln_new = n;
    mPolc = c;

    % frictions 
    mLambda_new(mPolaprime_new>p.Mina) = 0;
    mPolaprime_new(mPolaprime_new<=p.Mina) = p.Mina;

    %==========================
    % COMPUTE STATIONARY DISTRIBUTION
    %==========================

    % non-stochastic iterative histogram method
    iter2 = 1;
    err2 = 10;
    while err2 > tolDist

        mNewDist = zeros(size(mCurrDist));
        for iz = 1:p.Nz 
            for ia = 1:p.Na
                
                % for state (ia,iz) interpolate optimal savings on asset grid
                aprime = mPolaprime_new(ia,iz);
                nLow = sum(vGrida < aprime);
                nLow(nLow<=1) = 1;
                nLow(nLow>=p.Na) = p.Na - 1;
                nHigh = nLow + 1;
                wtLow = (vGrida(nHigh)-aprime)/(vGrida(nHigh)-vGrida(nLow));
                wtLow(wtLow>1) = 1;
                wtLow(wtLow<0) = 0;
                wtHigh = 1-wtLow;
                
                % get current mass over state (ia,iz) 
                mass = mCurrDist(ia,iz);
                
                % update mass cumutatively for all possible future iz states
                for izprime = 1:p.Nz
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
        iter2 = iter2+1;
    end
    
    %==========================
    % COMPUTE AGGREGATES
    %==========================

    % market clearing 
    margdist            = sum(mCurrDist,2);
    K_new               = vGrida * margdist;
    L_new               = sum(mgridz.*mPoln_new.*mCurrDist, 'all');
    avgmLambda          = sum(mLambda.*mCurrDist, 'all'); 
    avgmLambda_new      = sum(mLambda_new.*mCurrDist, 'all');
    avgmPolaprime       = sum(mPolaprime.*mCurrDist, 'all');
    avgmPolaprime_new   = sum(mPolaprime_new.*mCurrDist, 'all');
    avgmPoln            = sum(mPoln.*mCurrDist, 'all');
    avgmPoln_new        = sum(mPoln_new.*mCurrDist, 'all');

    %==========================
    % UPDATING
    %==========================
    
    % error 
    err1 = max(abs([...
        K_new               - K;...
        L_new               - L;...
        avgmLambda_new      - avgmLambda;...
        avgmPolaprime_new   - avgmPolaprime;...
        avgmPoln_new        - avgmPoln])); 

    % smooth updating 
    K           = wtOld1*K              + (1-wtOld1)*K_new;
    L           = wtOld2*L              + (1-wtOld2)*L_new;
    mPolaprime  = wtOld3*mPolaprime     + (1-wtOld3)*mPolaprime_new;
    mPoln       = wtOld4*mPoln          + (1-wtOld4)*mPoln_new;
    mLambda     = wtOld5*mLambda        + (1-wtOld5)*mLambda_new;

    %==========================
    % PROGRESS REPORTING
    %==========================

    timer = toc;
    if mod(iter1, dispiter) == 0
        fprintf('-----------------------------------------------------------------\n');
        fprintf('market clearing results: iteration %d in %.2fs\n', iter1, timer);
        fprintf('error: %.12f \n', err1);
        fprintf('K: %.6f    L: %.6f     r: %.6f    w: %.6f\n', K, L, r, w);
        fprintf('min n: %.3f    max n: %.3f   min c: %.3f    max c: %.3f\n',...
            min(min(mPoln)), max(max(mPoln)), min(min(mPolc)), max(max(mPolc)));
        fprintf('error lambda: %.8f    min lambda: %.3f    max lambda: %.3f\n',...
            abs(avgmLambda_new-avgmLambda), min(min(mLambda)), max(max(mLambda)));
        fprintf('Dist Converged in %d iterations\n', iter2);
        fprintf('\n');
        fprintf('err K: %.8f    err L: %.8f \n err aprime: %.8f    err n: %.8f\n',abs(K_new-K), abs(L_new-L),...
            max(abs(avgmPolaprime_new-avgmPolaprime)), max(abs(avgmPoln_new-avgmPoln)))
    end
    iter1 = iter1 + 1;
end


