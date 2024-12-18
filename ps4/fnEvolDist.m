%% Updating Joint Distribution Over (a,z) States %%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Update Joint-Distribution According to G_{t}(\Phi_{t}) = \Phi_{t+1}
%
%   Args:
%       params: model specific parameters 
%       mCurrDist: (matrix) current joint distribution over states (a,z)
%       mPolaprime: (matrix) current savings policy rule (interpolated)
%
%   Returns:
%       mNewDist: (matrix) updated joint distribution over states (a,z)
%
function [mNewDist] = fnEvolDist(params, opts_fnEvolDist)
    % unpack parameters
    vGrida2 = params.vGrida2;
    pNz = params.Nz;
    pNa2 = params.Na2;
    mPz = params.mPz;

    % unpack options 
    mCurrDist = opts_fnEvolDist.mCurrDist;
    mPolaprime = opts_fnEvolDist.mPolaprime;

    % reset new distribution (use Histogram method to incrementally update mass)
    mNewDist = zeros(size(mCurrDist));
    
    % update distribution using current distribution, savings policy rule and z-transition matrix 
    for iz = 1:pNz
        for ia = 1:pNa2
            aprime = mPolaprime(ia, iz);
            [L, H, wtL, wtH] = fnInterp1dGrid(aprime, vGrida2, pNa2);
            mass = mCurrDist(ia, iz);

            % updating distribution using \Phi_{t+1} = G_{t}(\Phi_{t}
            for iznext = 1:pNz
                mNewDist(L, iznext) = mNewDist(L, iznext) ...
                    + mPz(iz, iznext) * mass * wtL;

                mNewDist(H, iznext) = mNewDist(H, iznext) ...
                    + mPz(iz, iznext) * mass * wtH;
            end
        end
    end
end