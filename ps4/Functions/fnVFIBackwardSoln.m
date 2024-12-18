%% Value Function and Policy Rule Computaiton %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Single iteraiton solving the HH problem for given prices and future VF
%
%   Args:
%       params: fixed parameters defining the household problem
%       mVF: (matrix) future (next period) value function
%       aggK: (scalar) current aggregate capital allocation
%       aggL: (scalar) current aggregate labour supply allocation
%
%   Returns:
%       mVFnew: (matrix) current period value function (VF(t))
%       mPolaprime2: (matrix) current period savings policy rule (g_{a,t})
%       mPolc: (matrix) current period consumption policy rule (g_{c,t})
%
function [mVFnew, mPolaprime2, mPolc] = fnVFIBackwardSoln(params, opts_fnVFIBackwardSoln)
    % unpack options 
    mVF = opts_fnVFIBackwardSoln.VF;
    aggK = opts_fnVFIBackwardSoln.aggK;
    aggL = opts_fnVFIBackwardSoln.aggL;
    
    % unpack parameters 
    mPz = params.mPz;
    vGridz = params.vGridz;
    vGrida1 = params.vGrida1;
    vGrida2 = params.vGrida2;
    pNa1 = params.Na1;
    pNa2 = params.Na2;
    pNz = params.Nz;
    MinVala = params.minwealth;
    pAlpha = params.alpha;
    pTFP = params.tfp;
    pDelta = params.delta;
    pBeta = params.beta;

    % set placeholder variables 
    mVFnew = zeros(pNa1, pNz);
    mPolaprime1 = zeros(pNa1, pNz);
    mPolaprime2 = zeros(pNa2, pNz);
    mPolc = zeros(pNa1, pNz);

    % update prices given aggregate allocation
    r = pAlpha * pTFP * (aggK / aggL)^(pAlpha - 1) - pDelta;
    w = (1 - pAlpha) * pTFP * (aggK / aggL)^pAlpha;

    % compute VF and Policy Rules given next period VF and prices
    for iz = 1:pNz
        z = vGridz(iz);
        expVF = mVF * mPz'; % expected future value over all current states 
        expVF = expVF(:, iz); % expected future value conditional on current state iz
        minWealth = MinVala; % reset lower bound for current state
        for ia = 1:pNa1 
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
        end
    end

    % interpolate asset Policy Rule 
    for iz = 1:pNz
    mPolaprime2(:, iz) = interp1(vGrida1, mPolaprime1(:, iz), vGrida2, "linear", "extrap");
    end
end