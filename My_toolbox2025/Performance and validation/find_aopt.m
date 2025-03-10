function Aopt = find_aopt(RMSECV, pf, method)
% ----- Automatated method using punishing factor for choice of Aopt in PLSR models ------
% INPUTS:
%       RMSECV - the RMSECV output from from crossvalidation (for all #comps included)
%       pf     - Optional input, punishing factor
%       method - Optional input, 1 or 2 (1: Westad&Martens (default), 2:Alternative method*)
%                 *NB! method = 2 assumes that first elements in RMSECV is the residual 
%                      variance using zero components
% -------------------------------------------------------------------------
% OUTPUTS:
%       Aopt - scalar, the optimal number of PLS components to include in model
% -------------------------------------------------------------------------
% Option 1 method:
% Reference: F. Westad and H. Martens, 2000, NIR publications
% https://www.osapublishing.org/jnirs/abstract.cfm?uri=jnirs-8-2-117
%
% Option 2 method: 
% Alternative method
% -------------------------------------------------------------------------
% 
% Room for improvement: Input parsing
%
% Revision 12.01.22 : Westad and Martens generalized to not only look at
% previous number of components in the reduction if test.
% -------------------------------------------------------------------------

if nargin > 1
    PunishFactor = pf;
else
    PunishFactor = 0.03; % 3 percent
end

if nargin > 2 
    method = method;
else
    method = 1 ; 
end


if method == 1
    % Method 1: Westad&Martens -----------------------------------------------
    [~, Aopt] = min(RMSECV); % Find AOpt from minimum(RMSECV)for all no. components
    ReduceOneMore = 'yes';
    while strcmp(ReduceOneMore,'yes') 
        if Aopt == 1
            ReduceOneMore = 'no'; % Prevent reducing #LVs to negative values 
        elseif any(RMSECV(Aopt)*(1+PunishFactor) > RMSECV(1: Aopt-1 ))     
            candidates = find(RMSECV(Aopt)*(1+PunishFactor) > RMSECV(1: Aopt-1 ));
            Aopt = candidates(1); % Choose the lowest number of components amongst the candidates.
        else
            ReduceOneMore = 'no';
        end
    end
    
else
    % Method 2: Alternative method ----------------------------------------
    warning('Assumes that first elements in RMSECV input vector is the residual variance using zero components')
    A = length(RMSECV) -1 ;    
    crit = [0:A]*0.02*RMSECV(1) + RMSECV;
    [~,Aopt]= min(crit(2:end));
    % ---------------------------------------------------------------------
end

% 



end