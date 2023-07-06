function [dt] = ProcessRawDataNoSR(params, dt, path)
% Calculateand smooth raw velocities
dt.Vr = params.SampRateTreadmill*smooth(-0.5 * params.CalibRot * (dt.X1 + dt.X2), ...
    params.SmoothParameterBehavior/length(dt.X1), 'lowess');
dt.Vf = params.SampRateTreadmill*smooth(-0.5 * params.CalibForw * sqrt(2) * (dt.Y1 + dt.Y2), ...
    params.SmoothParameterBehavior/length(dt.X1), 'lowess');
dt.Vs = params.SampRateTreadmill*smooth(-0.5 * params.CalibForw * sqrt(2) * (dt.Y1 - dt.Y2), ...
    params.SmoothParameterBehavior/length(dt.X1), 'lowess');
% If closed loop align rotation to visual stimulus
if ~strcmp(params.pType, 'MovShort') && ~strcmp(params.pType, 'Dark1Min')
    dt.VisualStimulusCL = dt.Vr(dt.displayT);
else
    dt.displayT = zeros(size(dt.Vr));
end
% Crop membrane potential and photodiode signal to match the behavior data
dt.MembranePotential = dt.MembranePotential(params.TriggerTStart:params.TriggerTEnd);
dt.FlickerSignal = dt.FlickerSignal(params.TriggerTStart:params.TriggerTEnd);

% Calculate Protocol Transitions
tTransitionDAC = 0;
tTreadmill = dt.ID/params.SampRateTreadmill;
tTransitionTreadmill = 0;

for i = 1 : (length(params.ProtocolTransitionTime)-1)
    tTransSoft = find(dt.SoftwareTime >= params.ProtocolTransitionTime(i+1), 1);
    if(params.FlickerON == 1) % if using photodiode
        IndexFlickerPH = dt.ID(tTransSoft)/params.PhotodiodeFlickRate; % get the transition time in photdiode units
        if length(dt.PhotodiodeChange) > floor(IndexFlickerPH) % if transition time is smaller than length of trial
            tPreviousPHFlicker = dt.PhotodiodeChange(floor(IndexFlickerPH)); % find the flicker that preceeded the transition time 
            tTransitionDAC = vertcat(tTransitionDAC, tPreviousPHFlicker + ...
                floor((IndexFlickerPH-floor(IndexFlickerPH))*mean(params.AverageFlickTime))); % find the estimate transition time in DAC units
        else 
            tTransitionDAC = vertcat(tTransitionDAC, dt.PhotodiodeChange(end) + ...
                floor((IndexFlickerPH-floor(IndexFlickerPH))*mean(params.AverageFlickTime)));
        end
    else
        tTransitionDAC = vertcat(tTransitionDAC, tTreadmill(tTransSoft)*params.SampRateDAC);
    end
    tTransitionTreadmill = vertcat(tTransitionTreadmill, floor(params.SampRateTreadmill*tTransitionDAC(end)/params.SampRateDAC));
end

%if the protocol used is 'MovMultiple3' make the proper visual stimulations
if ~strcmp(params.pType, 'MovShort') && ~strcmp(params.pType, 'Dark1Min')
    if strcmp(params.pType, 'MovMultiple3')
        dt.VisualStimulus = zeros(size(dt.Vr));
        for k = 2 : length(tTransitionTreadmill)
            if mod(k,2) == 0
                dt.VisualStimulus((1+tTransitionTreadmill(k-1)):tTransitionTreadmill(k)) = 1;
            else
                dt.VisualStimulus((1+tTransitionTreadmill(k-1)):tTransitionTreadmill(k)) = 0;
            end
        end
    else
        dt.VisualStimulus = zeros(size(dt.Vr));
        for k = 2 : length(tTransitionTreadmill)
            if mod(k-1,4) == 2
                dt.VisualStimulus((1+tTransitionTreadmill(k-1)):tTransitionTreadmill(k)) = 1;
            elseif mod(k-1,4) == 0
                dt.VisualStimulus((1+tTransitionTreadmill(k-1)):tTransitionTreadmill(k)) = -1;
            else
                dt.VisualStimulus((1+tTransitionTreadmill(k-1)):tTransitionTreadmill(k)) = 0;
            end
        end
    end
end

% Detect activity events
[dt.actState, dt.ActBouts] = GetBouts(dt.Vr, dt.Vf, dt.Vs);

ddd = diff(params.ProtocolTransitionTime);
dt.Data = cell(length(tTransitionDAC)-1,1);
% Crop the time series into chunks corresponding to the subtrials
for k = 1 : length(dt.Data)
    dt.Data{k}.ID = dt.ID((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dt.Data{k}.TX1 = dt.X1((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dt.Data{k}.TX2 = dt.X2((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dt.Data{k}.TY1 = dt.Y1((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dt.Data{k}.TY2 = dt.Y2((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dt.Data{k}.Vr = dt.Vr((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dt.Data{k}.Display = dt.displayT((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    if strcmp(params.pType,'CircRDts')
        if mod(k,2) == 1
            % If natural gain
            dt.Data{k}.VisCL = dt.VisualStimulusCL((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
        else
            % If reverse gain
            dt.Data{k}.VisCL = -dt.VisualStimulusCL((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
        end
    end
    dt.Data{k}.Vf = dt.Vf((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dt.Data{k}.Vs = dt.Vs((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    [dt.Data{k}.actState, dt.Data{k}.ActBouts] = GetBouts(dt.Data{k}.Vr, dt.Data{k}.Vf, dt.Data{k}.Vs);
    dt.Data{k}.Frames = dt.CameraFrames((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dt.Data{k}.SoftTime = dt.SoftwareTime((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dt.Data{k}.FlickerSignalO = dt.FlickerSignal(((tTransitionDAC(k)+1):tTransitionDAC(k+1)));
    dt.Data{k}.MembranePotentialO = dt.MembranePotential(((tTransitionDAC(k)+1):tTransitionDAC(k+1)));
    
    % Resize Flicker and Membrane potential to same time base as behavior
    tx = linspace(1,length(dt.Data{k}.FlickerSignalO), length(dt.Data{k}.FlickerSignalO));
    dt.Data{k}.FlickerSignal = resample(dt.Data{k}.FlickerSignalO,tx,...
        length(dt.Data{k}.Vr)/length(dt.Data{k}.FlickerSignalO), 'spline');
    dt.Data{k}.MembranePotential = resample(dt.Data{k}.MembranePotentialO,tx,...
        length(dt.Data{k}.Vr)/length(dt.Data{k}.MembranePotentialO), 'spline');

    % Resample Signals to the same time base 
    % Use intervals of 125ms for resampling
    delta = 500;
    ndiv = ceil((ddd(k)*params.SampRateTreadmill)/delta);
    rdiv = ceil(length(dt.Data{k}.Vr)/ndiv);
    mdiv = ceil(length(dt.Data{k}.MembranePotential)/ndiv);

    dt.Data{k}.VrRes = [];
    dt.Data{k}.DispRes = [];
    dt.Data{k}.VisCLRes = [];
    dt.Data{k}.VfRes = [];
    dt.Data{k}.VsRes = [];
    dt.Data{k}.acStRes = [];
    dt.Data{k}.MPRes = [];
    dt.Data{k}.FlickRes = [];
    if rdiv ~= 0
        for i = 1 : ndiv
            if (1+(i-1)*rdiv)<=length(dt.Data{k}.Vr)
                if  i ~= ndiv && i*rdiv < length(dt.Data{k}.Vr)
                    dt.Data{k}.VrRes = vertcat(dt.Data{k}.VrRes, resample(dt.Data{k}.Vr((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    dt.Data{k}.DispRes = vertcat(dt.Data{k}.DispRes, resample(dt.Data{k}.Display((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    if strcmp(params.pType,'CircRDts')
                        dt.Data{k}.VisCLRes = vertcat(dt.Data{k}.VisCLRes, resample(dt.Data{k}.VisCL((1+(i-1)*rdiv):(i*rdiv)), ...
                            delta, rdiv, 0, 0));
                    end
                    dt.Data{k}.VfRes = vertcat(dt.Data{k}.VfRes, resample(dt.Data{k}.Vf((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    dt.Data{k}.VsRes = vertcat(dt.Data{k}.VsRes, resample(dt.Data{k}.Vs((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    dt.Data{k}.acStRes = vertcat(dt.Data{k}.acStRes, resample(dt.Data{k}.actState((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    dt.Data{k}.MPRes = vertcat(dt.Data{k}.MPRes, resample(dt.Data{k}.MembranePotential((1+(i-1)*mdiv):(i*mdiv)), ...
                        delta, mdiv, 0, 0));
                    dt.Data{k}.FlickRes = vertcat(dt.Data{k}.FlickRes, resample(dt.Data{k}.FlickerSignal((1+(i-1)*mdiv):(i*mdiv)), ...
                        delta, mdiv, 0, 0));
                else
                    auxrsp = resample(dt.Data{k}.Vr((1+(i-1)*rdiv):end), delta, length(dt.Data{k}.Vr((1+(i-1)*rdiv):end)), 0, 0);
                    dt.Data{k}.VrRes = vertcat(dt.Data{k}.VrRes, auxrsp(1:end));
                    auxrsp = resample(dt.Data{k}.Display((1+(i-1)*rdiv):end), delta, length(dt.Data{k}.Display((1+(i-1)*rdiv):end)), 0, 0);
                    dt.Data{k}.DispRes = vertcat(dt.Data{k}.DispRes, auxrsp(1:end));
                    if strcmp(params.pType,'CircRDts')
                        auxrsp = resample(dt.Data{k}.VisCL((1+(i-1)*rdiv):end), delta, length(dt.Data{k}.VisCL((1+(i-1)*rdiv):end)), 0, 0);
                        dt.Data{k}.VisCLRes = vertcat(dt.Data{k}.VisCLRes, auxrsp(1:end));
                    end
                    auxrsp = resample(dt.Data{k}.Vf((1+(i-1)*rdiv):end), delta, length(dt.Data{k}.Vf((1+(i-1)*rdiv):end)), 0, 0);
                    dt.Data{k}.VfRes = vertcat(dt.Data{k}.VfRes, auxrsp(1:end));
                    auxrsp = resample(dt.Data{k}.Vs((1+(i-1)*rdiv):end), delta, length(dt.Data{k}.Vs((1+(i-1)*rdiv):end)), 0, 0);
                    dt.Data{k}.VsRes = vertcat(dt.Data{k}.VsRes, auxrsp(1:end));
                    auxrsp = resample(dt.Data{k}.actState((1+(i-1)*rdiv):end), delta, length(dt.Data{k}.actState((1+(i-1)*rdiv):end)), 0, 0);
                    dt.Data{k}.acStRes = vertcat(dt.Data{k}.acStRes, auxrsp(1:end));
                    auxrsp = resample(dt.Data{k}.MembranePotential((1+(i-1)*mdiv):end), delta, length(dt.Data{k}.MembranePotential((1+(i-1)*mdiv):end)), 0, 0);
                    dt.Data{k}.MPRes = vertcat(dt.Data{k}.MPRes,  auxrsp(1:end));
                    auxrsp = resample(dt.Data{k}.FlickerSignal((1+(i-1)*mdiv):end), delta, length(dt.Data{k}.FlickerSignal((1+(i-1)*mdiv):end)), 0, 0);
                    dt.Data{k}.FlickRes = vertcat(dt.Data{k}.FlickRes,  auxrsp(1:end));
                end
            end
        end
    end
end
% Save preprocess data
dt.params = params;
disp(['Transitions Segmented: ' num2str(length(tTransitionDAC)-1)]);
save([path '\PreProcessedDataNoSR.mat'], 'dt','-v7.3')
disp('PreProcessedDataNoSR.mat file created')

end