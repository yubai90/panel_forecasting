function [logPL,CRPS] = densfore(Y_forecast_sim,Y_f,h,NG)
% density forecast, used for parallel loop
logPL=zeros(NG,h);
CRPS=zeros(NG,h);
for i=1:NG
    logPL(i,:) = get_logscoreH(squeeze(Y_forecast_sim(:,i,:)), Y_f(:,i), 1:h, 'plain');
    CRPS(i,:) = get_CRPSH(squeeze(Y_forecast_sim(:,i,:)), Y_f(:,i), 1:h);
end

end

