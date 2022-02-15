function vecLogScores = get_logscoreH(matPaths, vecData, vecHorizons, strType)
% -------------------------------------------------------------------------
% get_logscore:
%
% Calculates log-scores for a set of predicted values and actual
% observations on a specific variable
% 
% INPUTS:   - matPaths: (# draws)*(# periods) matrix of simulated values for
%               the variable.
%           - vecData: (# periods)*1 vector of observations on the variable
%           - vecHorizons: vector defining the horizons at which the
%               forecasts are evaluated
%           - strType: 'plain' to score period-specific forecasts, 'cumul' to
%               score cumulative forecasts. Choice depends on the variable:
%               use 'cum' for changes/growth rates.
%
% OUTPUTS:  - vecLogScores: 1*... vector of log-scores. The value
%               is NaN for the densities for which we do not have
%               observations (ie the last ones).
%
% P Alessandri, Jan 2013
% Y BAI, compute log score by quadratic approximations
% -------------------------------------------------------------------------


% Initialise vector of log-scores:
vecLogScores = NaN(1,length(vecHorizons));

% Loop over horizons to get log-score for each predictive density

for i = 1:length(vecHorizons)
    
    ii=vecHorizons(i);
    
    % Define paths and observations as plain or cumulative:
    switch strType
        case 'plain'
            paths = matPaths(:,ii);
            obs   = vecData(ii);
        case 'cumul'
            paths = sum(matPaths(:,1:ii), 2);
            obs   = sum(vecData(1:ii));
    end
    
  %  vecLogScores(i)= -0.5*(log(2*pi)+log(var(paths))+(obs-mean(paths))^2/var(paths));
    
    vecLogScores(i) = log( ksdensity(paths, obs)+1e-9 );

end

end