%% Create Joint-Transition Matrix
function [mJointTrans] = fnTransMat(pNa, pNz, mPolaprime, vGrida, mPz)
% 
% Compute Joint-Transition Matrix Over (a,z) States 
%
% Args:
%   pNa: number of states for a (assets)
%   pNz: Number of states for z (productivity)
%   mPolaprime: (interpolated) policy rule for optimal savings choice 
%   vGrida: asset grid 
%   mPz: probability transition matrix for productivity
%
mJointTrans = zeros(pNa*pNz, pNa*pNz);
for ia = 1:pNa
    for iz = 1:pNz
        iCurrState = (iz - 1)*pNa + ia; % corresponding i-th row for state (ia, iz)
        aprime = mPolaprime(ia, iz); % optimal savings given current state (ia, iz)
        [LB, UB, wtLB, wtUB] = fnInterp1dGrid(aprime, vGrida, pNa); % interpolate aprime on asset grid
        
        for izfuture = 1:pNz
            iFutureStateLB = (izfuture - 1)*pNa + LB; % index for LB - future state 
            iFutureStateUB = (izfuture - 1)*pNa + UB; % index for UB - future state 
            mJointTrans(iCurrState, iFutureStateLB) = 1 * mPz(iz, izfuture) * wtLB; % prob. transition to LB
            mJointTrans(iCurrState, iFutureStateUB) = 1 * mPz(iz, izfuture) * wtUB; % prob. transition to UB
        end
    end
end