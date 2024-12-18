%% Solving Aiyagari (1994) Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Solving SRCE w/ Howard Improvement in VFI and Interpolated Policy Rule
%
%   Args:
%       alpha: (scalar) capital share of income 
%       beta: (scalar) discount factor (annual)
%       delta: (scalar) depreciation rate of capital
%       tfp: (scalar) TFP (A)
%       mu: (scalar) constant term in AR(1) process 
%       rho: (scalar) coefficient of persistence in AR(1) process
%       sigma: (scalar) scalar term multiplying white noise term in AR(1) process: 
%           y(t) = mu + rho*y(t-1) + sigma*e(t), e(t) ~ iid N(0,1)
%       Nz: (scalar) number of states to discretize AR(1) process (using Tauchen)
%       minwealth: (scalar) lower bound on savings (minimum wealth constraint)
%       initk: (scalar) initial aggregate K (guess) to start GE loop
%       progreport: (boolean) true - report loop progress. false - don't report loop progress
%
%   Returns:
%       mPolc: (matrix) consumption policy rule
%       mPolaprime2: (matrix) savings policy rule (interpolated)
%       mVF: (matrix) value function 
%       mCurrDist: (matrix) joint-distribution over states (a,z) 
%       aggK: (scalar) aggregate capital supply/demand
%       aggL: (scalar) aggregate labour supply
%       w: (scalar) equilibrium wage rate 
%       r: (scalar) equilibrium real rate 
%       vGrida1: (vector) sparse asset grid 
%       vGrida2: (vector) dense asset grid
%       vGridz: (vector) productivity state grid
%       mPz: (matrix) productivity transition matrix
%       vPiz: (vector) invariant distribution of markov process under mPz
%
function [mPolc, mPolaprime2, mVF, mCurrDist, aggK, aggL, w, r, vGrida1, vGrida2, vGridz, mPz, vPiz] ...
        = fnSolveSRCE(params)
    %% Set Parameters 
    % Aiyagari (1994) 
    pAlpha = params.alpha;
    pBeta = params.beta;
    pDelta = params.delta;
    pTFP = params.tfp;
    
    % income process 
    pMu = params.mu;
    pRho = params.rho;
    pSigma = params.sigma; %0.2 * sqrt(1 - pRho^2);
    pNz = params.Nz;
    
    % wealth grid 
    pNa1 = 50;
    pNa2 = 100;
    pMinGrida = params.minwealth;
    pMaxGrida = 150;
    pCurve = 7;
    
    % loop parameters 
    pWtOldK = 0.9900;
    pTolGE = 1e-8;
    pTolVF = 1e-8;
    pTolDist = 1e-8;
    pMaxGEiter = 2000;
    pMaxVFiter = 2000;
    convergedGE = false;
    progreport = params.progreport;
    %convergedVF = false;
    %% Define Grids 
    % wealth - coarsed grid (Maliar, Maliar and Valli, 2010)
    x1 = linspace(0, 0.5, pNa1);
    x2 = linspace(0, 0.5, pNa2);
    y1 = x1.^pCurve / max(x1.^pCurve);
    y2 = x2.^pCurve / max(x2.^pCurve);
    vGrida1 = pMinGrida + (pMaxGrida - pMinGrida) * y1;
    vGrida2 = pMinGrida + (pMaxGrida - pMinGrida) * y2;
    
    % productivity 
    [vZ, mPz] = fnTauchen(pRho, pSigma, pMu, pNz);
    vPiz = fnStationaryDist(mPz);
    vGridz = exp(vZ);
    %% Initialize 
    % initial guess (aggregate(s)) 
    aggKinit = params.initk;
    mVF = repmat(0.01 .* vGrida1', 1, pNz);
    mCurrDist = ones(pNa2, pNz);
    mCurrDist = mCurrDist ./ (pNa2 * pNz);
    
    % placeholders 
    mVFnew = zeros(pNa1, pNz);
    mPolaprime1 = zeros(pNa1, pNz);
    mPolc = zeros(pNa1, pNz);
    mPolaprime2 = zeros(size(mCurrDist));
    %aggKnew = 0;
    %mNewDist = zeros(pNa2, pNz);
    %% Solve Model
    aggK = aggKinit;
    aggL = vGridz' * vPiz; 
    iterGE = 1;
    errGE = 20;
    tic; 
    while errGE > pTolGE && iterGE <= pMaxGEiter
        r = pAlpha * pTFP * (aggK / aggL)^(pAlpha - 1) - pDelta;
        w = (1 - pAlpha) * pTFP * (aggK / aggL)^pAlpha;
        % --------------------------------------------------- %
        % INNER LOOP - VFI + HOWARD IMPROVEMENT
        % --------------------------------------------------- %
        iterVF = 1;
        errVF = 20;
        while errVF > pTolVF && iterVF <= pMaxVFiter
                    
            for iz = 1:pNz
    
                z = vGridz(iz);
                expVF = mVF * mPz'; % expected future value over all current states 
                expVF = expVF(:, iz); % expected future value conditional on current state iz
                minWealth = pMinGrida; % reset lower bound for current state 
        
                for ia = 1:pNa1
                    
                    % ----- VFI w/ Interpolated Policy Function ----- %        
                    if iterVF <= 30 || mod(iterVF, 20) == 0
                        a = vGrida1(ia);
                        budget = w*z + (1+r)*a; 
                        lb = minWealth;
                        ub = budget;
                        aprime = fnOptaprime({pBeta, budget, vGrida1, pNa1, expVF}, lb, ub);
                        c = budget - aprime;
        
                        [LB, UB, wtLB, wtUB] = fnInterp1dGrid(aprime, vGrida1, pNa1);
                        value = wtLB*expVF(LB) + wtUB*expVF(UB); % interpolate E[V]
                        
                        % updating
                        mVFnew(ia, iz) = log(c) + pBeta * value;
                        mPolaprime1(ia, iz) = aprime;
                        mPolc(ia, iz) = c;
                        minWealth = aprime; % exploit monotonicity of policy rule
                    
                    % ----- VFI w/ Accelerated Howard Improvement ----- % 
                    else    
                        aprime = mPolaprime1(ia, iz);
                        c = mPolc(ia, iz);
                        [LB, UB, wtLB, wtUB] = fnInterp1dGrid(aprime, vGrida1, pNa1);
                        value = wtLB*expVF(LB) + wtUB*expVF(UB); 
                        
                        % updating
                        mVFnew(ia, iz) = log(c) + pBeta * value;
                    end
                end
            end    
            errVF = max(max(abs(mVFnew - mVF)));
            mVF = mVFnew;    
            iterVF = iterVF + 1;
        end    
        % --------------------------------------------------- %
        %  COMPUTE STATIONARY DISTRIBUTION 
        % --------------------------------------------------- %    
        for iz = 1:pNz
            mPolaprime2(:, iz) = interp1(vGrida1, mPolaprime1(:, iz), vGrida2, "linear", "extrap");
        end
        errDist = 20;
        iterDist = 1;
        while errDist > pTolDist
    
            mNewDist = zeros(size(mCurrDist)); % reset distribution to incremeentally increase mass over states (ia,iz)
    
            for iz = 1:pNz
                for ia = 1:pNa2    
                    aprime = mPolaprime2(ia, iz);
                    [L, H, wtL, wtH] = fnInterp1dGrid(aprime, vGrida2, pNa2);
                    mass = mCurrDist(ia, iz);
                    
                    % updating
                    for iznext = 1:pNz
                        mNewDist(L, iznext) = mNewDist(L, iznext) ...
                            + mPz(iz, iznext) * mass * wtL;
    
                        mNewDist(H, iznext) = mNewDist(H, iznext) ...
                            + mPz(iz, iznext) * mass * wtH;
                    end
                end
            end
            errDist = max(max(abs(mNewDist - mCurrDist)));
            mCurrDist = mNewDist;
            iterDist = iterDist + 1;
        end
        % --------------------------------------------------- %
        %  UPDATE AGGREGATES USING MARKET CLEARING 
        % --------------------------------------------------- %
        vMarginalDista = sum(mCurrDist, 2); % (marginal) density over asset-grid
        Kmcc = vGrida2 * vMarginalDista; % capital market clearing condition
        errGE = abs(Kmcc - aggK);
    
        % updating
        aggKnew = pWtOldK*aggK + (1-pWtOldK)*Kmcc; % smooth updating     
        aggK = aggKnew; 
        timer = toc;
        % --------------------------------------------------- %
        %  REPORT PROGRESS 
        % --------------------------------------------------- %
        if mod(iterGE, 20) == 0 && progreport
            fprintf('Iteration: %d in %.3f mins. Error: %.8f\n', iterGE, timer./60, errGE);
            fprintf('AggK: %.4f\n', aggK);
            fprintf('Ea(r): %.4f\n', Kmcc)
            fprintf('r: %.4f, w: %.4f\n', r, w);
            fprintf('----------------------------------------------------\n')
        end
        iterGE = iterGE + 1;
    end
    % Convergence Report 
    if errGE <= pTolGE
        convergedGE = true;
    end
    fprintf('GE Converged: %s\n', mat2str(convergedGE));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%