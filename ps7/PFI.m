%% Housekeeping 
close all;clc;
clear variables;
% addpath('./Figures')
addpath('./Functions')

%% MODEL FUNDAMENTALS  
% model parameters 
p.Alpha     = 0.36;
p.Beta      = 0.99;
p.Delta     = 0.025;
p.Theta     = 0.50;
p.Frisch    = 1.00;
p.Eta       = 7.60;

% other parameters
p.Rhoz      = 0.90;
p.Sigmaz    = 0.05;
p.Nz        = 7;
p.Mina      = 1e-20;
p.Maxa      = 300;
p.Na        = 100;
p.Curve     = 7;

% unpack params 
pbeta=p.Beta;palpha=p.Alpha;pdelta=p.Delta;pfrisch=p.Frisch;peta=p.Eta;

%% GRIDS AND RANDOM PROCESSES 
% asset grid 
x = linspace(0,0.5,p.Na);
y = x.^p.Curve / max(x.^p.Curve);
vGrida = p.Mina + (p.Maxa - p.Mina).*y;

% idiosyncratic productivity grid
[vGridz, mPz] = fnTauchen(p.Rhoz,p.Sigmaz,0.0,p.Nz);
vGridz = exp(vGridz);

%% SRCE SETUP 
% auxillary objects 
mgrida = repmat(vGrida', 1, p.Nz); 
mgridz = repmat(vGridz', p.Na, 1);

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
L           = 0.33;
mPolaprime  = zeros(size(mPolc)); 
mPoln       = ones(size(mPolc)).*L;
mLambda     = zeros(size(mPolc));
mCurrDist   = ones(p.Na,p.Nz)/(p.Na*p.Nz);

% loop parameters 
tol_ge = 1e-8;
tol_pfi = 1e-8;
tol_dist = 1e-8;
wtOld1 = 0.99500;
wtOld2 = 0.99500;
wtOld3 = 0.90000;
wtOld4 = 0.90000;
wtOld5 = 0.90000;

%=========================
% GE LOOP
%=========================

iter1 = 1;
err1 = 10;
idisp_outer = 5;
idisp_inner = 2500;
tic;
while err1 > tol_ge 

    %=========================
    % UPDATE PRICES
    %=========================

    r = palpha * A * (K/L)^(palpha-1) - pdelta;
    w = (1-palpha) * A * (K/L)^palpha;
    mmu = r+pdelta;

    %=========================
    % INNER LOOP - PFI 
    %%=========================

    iter2 = 1;
    err2 = 10;
    timer_temp = toc;
    while err2 > tol_pfi 

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
            
            % cumutitively compute expectation term 
            cprime = wprime.*zprime.*nprime + (1+rprime).*mPolaprime - aprimeprime;
            cprime(cprime < 1e-10) = 1e-10;
            muprime = 1./cprime;
            mExp = mExp + repmat(mPz(:,izprime)', p.Na, 1).*(1+rprime).*muprime;
        end
    
        %==========================
        % OPTIMAL POLICY GIVEN BELIEFS
        %==========================
        
        % intratemporal choice variables 
        mExp = pbeta*mExp;
        c = 1./(mExp+mLambda);
        c(c < 1e-10) = 1e-10;
        mPoln_new = ((w*mgridz)./(peta*c)).^pfrisch;
        mPoln_new(mPoln_new > 1) = 1;
        mPoln_new(mPoln_new < 0) = 0;
        
        % savings rule and lagrangian constraint 
        mLambda_new = 1./((1+r)*mgrida + w*mgridz.*mPoln - mPolaprime) - mExp;
        mPolaprime_new = (1+r)*mgrida + w*mgridz.*mPoln - c;
       
        % frictions 
        mLambda_new(mPolaprime_new>p.Mina) = 0;
        mPolaprime_new(mPolaprime_new<=p.Mina) = p.Mina;

        % compute error 
        err2 = abs(sum(mPolaprime_new.*mCurrDist, 'all') - sum(mPolaprime.*mCurrDist, 'all'));
        %err2 = max(abs(mPolaprime_new-mPolaprime), [], "all");
        err_n = max(abs(mPoln_new-mPoln), [], "all");
        err_lambda = max(abs(mLambda_new-mLambda), [], "all");

        % smooth updating 
        mPolaprime  = wtOld3.*mPolaprime    + (1-wtOld3).*mPolaprime_new;
        mPoln       = wtOld4.*mPoln         + (1-wtOld4).*mPoln_new;
        mLambda     = wtOld5.*mLambda       + (1-wtOld5).*mLambda_new;
        mPolc = c;

        % track progress
        % if mod(iter2, idisp_inner) == 0
        %     fprintf(' \n')
        %     fprintf('PFI LOOP \n')
        %     fprintf('iteration %d   error: %.15f \n', iter2, err2)
        %     fprintf('err_n: %.15f   err_lambda: %.15f \n', err_n, err_lambda)
        %     fprintf(' \n')
        % end
        iter2 = iter2+1;
    end
    timer1 = toc-timer_temp;

    %==========================
    % COMPUTE STATIONARY DISTRIBUTION
    %==========================

    % non-stochastic iterative histogram method
    iter3 = 1;
    err3 = 10;
    timer_temp = toc;
    while err3 > tol_dist

        mNewDist = zeros(size(mCurrDist));
        for iz = 1:p.Nz 
            for ia = 1:p.Na
                
                % for state (ia,iz) interpolate optimal savings on asset grid
                aprime = mPolaprime(ia,iz);
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
        err3 = max(abs(mNewDist-mCurrDist),[],'all');
        mCurrDist = mNewDist;
        iter3 = iter3+1;
    end
    timer2 = toc-timer_temp;
    
    %==========================
    % COMPUTE AGGREGATES
    %==========================

    % market clearing 
    margdist            = sum(mCurrDist,2);
    K_new               = vGrida * margdist;
    L_new               = sum(mgridz.*mPoln.*mCurrDist, 'all');
    avgmLambda          = sum(mLambda.*mCurrDist, 'all'); 
    avgmPolaprime       = sum(mPolaprime.*mCurrDist, 'all');
    avgmPoln            = sum(mPoln.*mCurrDist, 'all');

    %==========================
    % UPDATING
    %==========================
    
    % error 
    err1 = mean(abs([...
        K_new - K;...
        L_new - L]),...
        'all'); 

    % smooth updating 
    K = wtOld1*K + (1-wtOld1)*K_new;
    L = wtOld2*L + (1-wtOld2)*L_new;
   
    %==========================
    % PROGRESS REPORTING
    %==========================

    timer3 = toc;
    if mod(iter1, idisp_outer) == 0 || iter1 == 1
        fprintf('-----------------------------------------------------------------\n');
        fprintf('market clearing results: \n') 
        fprintf('iteration %d in %.2fs    error: %.15f \n', iter1, timer3, err1);
        fprintf('K: %.12f \n', K)
        fprintf('L: %.12f \n', L)
        fprintf('r: %.12f \n', r)
        fprintf('w: %.12f \n', w)
        fprintf('min n: %.4f    max n: %.4f   min c: %.4f    max c: %.4f\n', min(min(mPoln)), max(max(mPoln)), min(min(mPolc)), max(max(mPolc)))
        fprintf('min_lambda: %.6f    max_lambda: %.6f \n', min(min(mLambda)), max(max(mLambda)))
        fprintf('Dist Converged in %d iterations in %.2fs \n', iter3, timer2)
        fprintf('PFI converged in %d iterations in %.2fs \n', iter2, timer1)
    end
    
    
    % if iter1 >= 500
    %     % wtOld1 = 0.99900;
    %     % wtOld2 = 0.99900;
    %     wtOld3 = 0.99500;
    %     wtOld4 = 0.99500;
    %     wtOld5 = 0.99500;
    % end
    % elseif iter1 >= 5000 
    %     wtOld1 = 0.99950;
    %     wtOld2 = 0.99950;
    %     wtOld3 = 0.99500;
    %     wtOld4 = 0.99500;
    %     wtOld5 = 0.99500;
    % end
    iter1 = iter1 + 1;
end


