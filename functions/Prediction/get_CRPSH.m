function vecCRPS = get_CRPSH(matPaths, vecData, vecHorizons)
% -------------------------------------------------------------------------
% Compute continuous ranked probability scor 
% -------------------------------------------------------------------------


% Initialise vector of log-scores:
vecCRPS = NaN(1,length(vecHorizons));

% Loop over horizons to get log-score for each predictive density

for i = 1:length(vecHorizons)
    
    ii=vecHorizons(i);
    
    % Define paths and observations 
    paths = matPaths(:,ii);
    obs   = vecData(ii);
        
    vecCRPS(i) = crps(paths,obs);

end

end