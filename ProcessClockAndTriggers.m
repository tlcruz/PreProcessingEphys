function [params, dt] = ProcessClockAndTriggers(params, dt, VStimON)
% Clock and Trigger Signals, start looking only 1000 samples in the future
[VON, triggTBdacON] = max(diff(dt.TreadmillTrigger(1000:end)));   % find protocol start in ADC data
if abs(VON)<1
    disp('No trigger ON detected')
    triggTBdacON = 1;
    params.TrigONdetect = 'No trigger ON detected';
end
[VOFF, triggTBdacOFF] = min(diff(dt.TreadmillTrigger(1000:end)));  % find protocol end in ADC data
if abs(VOFF)<1
    disp('No trigger OFF detected')
    triggTBdacOFF = length(dt.TreadmillTrigger(1000:end))-10;
    params.TrigOFFdetect = 'No trigger OFF detected';
end
% Pass triggers to 1000 samples in the past
triggTBdacON = triggTBdacON + 1000;
triggTBdacOFF = triggTBdacOFF + 1000;

% Check if treadmill data was lost
dclockTB = diff(dt.TClock);
dclockTB(dclockTB==-254) = 1;
isloss = length(find(dclockTB~=1));
if(isloss == 0)
    disp('No Treadmill Data Lost');
    params.TreadmillLostData = 'No Treadmill Data Lost';
else
    disp(['Treadmill Data Lost ' num2str(isloss) ' times.']);
    params.TreadmillLostData = ['Treadmill Data Lost ' num2str(isloss) ' times.'];
end

% Display durations of all the different clocks
disp('---------------------------------------')
disp(['Set: ' num2str(params.ExperimentTime) ' s' ])
disp(['Treadmill: ' num2str(length(dt.TClock)/params.SampRateTreadmill) ' s' ])
disp(['DAC: ' num2str((triggTBdacOFF-triggTBdacON)/params.SampRateDAC) ' s' ])
disp(['Software: ' num2str(dt.SoftwareTime(end)-dt.SoftwareTime(1)) ' s' ])
disp('---------------------------------------')

auxFrT = length(dt.TClock)/params.SampRateTreadmill;
auxFrDAC = (triggTBdacOFF-triggTBdacON)/params.SampRateDAC;
params.SampRateDAC = params.SampRateDAC-params.SampRateDAC*(auxFrT-auxFrDAC)/auxFrDAC;
params.TriggerTStart = triggTBdacON;
params.TriggerTEnd = min(triggTBdacOFF,length(dt.FlickerSignal));

if VStimON == 1
    % If the flicker object is on extract the oscillation in photodiode
    % signal
    phFlick = zeros(length(dt.FlickerSignal),1);
    auxPh = smooth(dt.FlickerSignal, 25/length(dt.FlickerSignal), 'lowess');
    phCents = -1:0.005:1;
    [~, loc, ~, p] = findpeaks(-hist(auxPh, phCents));
    [~, I] = max(p);
    thr = phCents(loc(I));
    if thr  < 0.1 || thr  > 0.7
        thr = 0.4;
    end
    phFlick(auxPh > thr) = 1;
    % Measure the change on the photodiode signal
    phFlick = phFlick(params.TriggerTStart:params.TriggerTEnd);
    phChange = find(abs(diff(phFlick)) > 0);
    photodiodeFrameRate = phChange(2:end)-phChange(1:end-1);
    if (phChange(1) - mean(photodiodeFrameRate)) > 0
        phChange = vertcat(phChange(1)- mean(photodiodeFrameRate), phChange);
    end
    
    disp(['Target Flicker Time: ' num2str(params.PhotodiodeFlickRate/params.SampRateTreadmill) ' s'])
    disp(['Average Measured Flicker Time: ' num2str(mean(photodiodeFrameRate)/params.SampRateDAC) ' s'])
    disp(['Noise Flicker Time: ' num2str(100*std(photodiodeFrameRate)/mean(photodiodeFrameRate)) ' %'])
    disp('----------------------------------------')
    
    dt.PhotodiodeChange = phChange;
    params.FSQuality = 1;
    params.AverageFlickTime = mean(photodiodeFrameRate)/params.SampRateDAC;
    % Compare the flicker time with the preset and send warning if it does
    % there is a mismatch
    if abs(params.PhotodiodeFlickRate/params.SampRateTreadmill - mean(photodiodeFrameRate)...
            /params.SampRateDAC)>std(photodiodeFrameRate)/params.SampRateTreadmill
        disp('----------------------------------------')
        disp('WARNING - SOMETHING WRONG WITH FLICKER SIGNAL')
        disp('----------------------------------------')
        disp('----------------------------------------')
        params.FSQuality = 0;
        params.AverageFlickTime = params.PhotodiodeFlickRate/params.SampRateTreadmill;
        auxx = params.PhotodiodeFlickRate/params.SampRateTreadmill;
        avD = 0.04*params.SampRateDAC;
        dt.PhotodiodeChange = (avD+auxx*params.SampRateDAC:auxx*params.SampRateDAC:(auxx*params.SampRateDAC+params.SampRateDAC*params.ExperimentTime))';
    end
    
    % Compare treadmill clock with photodiode flicker
    tTimeFlick = (0:params.PhotodiodeFlickRate:length(dt.TClock))/params.SampRateTreadmill;
    tTimeFlick = tTimeFlick(1:min(length(tTimeFlick),length(phChange)))';
    phChangeTime = phChange(1:min(length(tTimeFlick),length(phChange)))/...
        params.SampRateDAC;
    % Measure lags in visual display
    delays = phChangeTime - tTimeFlick;
    disp(['Lag: ' num2str(1000*mean(delays)) ' +/- ' num2str(1000*std(delays)) ' ms  (FrameTime = ' num2str(1000/60) ') ms'])
    disp('----------------------------------------')
    
    % Align visual display given the instantaneous lags
    phChangeTMdisplay = phChangeTime*params.SampRateTreadmill;
    tTFlick = tTimeFlick*params.SampRateTreadmill;
    displayT = ones(size(dt.TClock));
    displayT(floor(phChangeTMdisplay)) = floor(tTFlick);
    inds = find(displayT>1);
    for k = 1 : length(inds)+1
        if k == 1
            lengI = inds(k);
            inp = round(linspace(displayT(1), displayT(inds(k)), lengI));
            displayT(1:inds(k)) = inp;
        elseif k > 1 && k < length(inds)+1
            lengI = inds(k)-inds(k-1);
            inp = floor(linspace(displayT(inds(k-1)), displayT(inds(k)), lengI));
            displayT((1+inds(k-1)):inds(k)) = inp;
        else
            lengI = length(displayT)-inds(end);
            inp = floor(linspace(displayT(inds(end)), length(displayT), lengI));
            displayT((1+inds(end)):length(displayT)) = inp;
        end
    end
    dt.displayT = displayT;
else
    dt.PhotodiodeChange = zeros(size(dt.FlickerSignal));
end
end













