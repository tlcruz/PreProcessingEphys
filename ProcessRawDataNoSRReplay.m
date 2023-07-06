function [dtR] = ProcessRawDataNoSRReplay(params, dtR, dtCL, path)
% Calculateand smooth raw velocities
dtR.Vr = params.SampRateTreadmill*smooth(-0.5 * params.CalibRot * (dtR.X1 + dtR.X2), ...
    params.SmoothParameterBehavior/length(dtR.X1), 'lowess');
dtR.Vf = params.SampRateTreadmill*smooth(-0.5 * params.CalibForw * sqrt(2) * (dtR.Y1 + dtR.Y2), ...
    params.SmoothParameterBehavior/length(dtR.X1), 'lowess');
dtR.Vs = params.SampRateTreadmill*smooth(-0.5 * params.CalibForw * sqrt(2) * (dtR.Y1 - dtR.Y2), ...
    params.SmoothParameterBehavior/length(dtR.X1), 'lowess');

% Align rotation to visual stimulus
dtR.VisualStimulusRep = nan*ones(size(dtR.displayT));
dtR.displayT(dtR.displayT > length(dtCL.Vr)) = length(dtCL.Vr);
if length(dtR.displayT) < length(dtCL.Vr)
    dtR.VisualStimulusRep(1:length(dtR.displayT)) = dtCL.Vr(dtR.displayT); 
else
    dtR.VisualStimulusRep(1:length(dtCL.Vr)) = dtCL.Vr(dtR.displayT(1:length(dtCL.Vr)));
end
% Crop membrane potential and photodiode signal to match the behavior data
dtR.MembranePotential = dtR.MembranePotential(params.TriggerTStart:params.TriggerTEnd);
dtR.FlickerSignal = dtR.FlickerSignal(params.TriggerTStart:params.TriggerTEnd);

% Calculate Protocol Transitions
tTransitionDAC = 0;
tTreadmill = dtR.ID/params.SampRateTreadmill;
tTransitionTreadmill = 0;
for i = 1 : (length(params.ProtocolTransitionTime)-1)
    tTransSoft = find(dtR.SoftwareTime >= params.ProtocolTransitionTime(i+1), 1);
    if(params.FlickerON == 1) % if using photodiode
        IndexFlickerPH = dtR.ID(tTransSoft)/params.PhotodiodeFlickRate; % get the transition time in photdiode units
        if length(dtR.PhotodiodeChange) > floor(IndexFlickerPH) % if transition time is smaller than length of trial
            tPreviousPHFlicker = dtR.PhotodiodeChange(floor(IndexFlickerPH)); % find the flicker that preceeded the transition time 
            tTransitionDAC = vertcat(tTransitionDAC, tPreviousPHFlicker + ...
                floor((IndexFlickerPH-floor(IndexFlickerPH))*mean(params.AverageFlickTime))); % find the estimate transition time in DAC units
        else 
            tTransitionDAC = vertcat(tTransitionDAC, dtR.PhotodiodeChange(end) + ...
                floor((IndexFlickerPH-floor(IndexFlickerPH))*mean(params.AverageFlickTime)));
        end
    else
        tTransitionDAC = vertcat(tTransitionDAC, tTreadmill(tTransSoft)*params.SampRateDAC);
    end
    tTransitionTreadmill = vertcat(tTransitionTreadmill, floor(params.SampRateTreadmill*tTransitionDAC(end)/params.SampRateDAC));
end
% Detect activity events
[dtR.actState, dtR.ActBouts] = GetBouts(dtR.Vr, dtR.Vf, dtR.Vs);
[dtR.visState, dtR.VisBouts] = GetBouts(dtR.VisualStimulusRep, zeros(size(dtR.Vf)), zeros(size(dtR.Vs)));

ddd = diff(params.ProtocolTransitionTime);
dtR.Data = cell(length(tTransitionDAC)-1,1);
% Crop the time series into chunks corresponding to the subtrials
for k = 1 : length(dtR.Data)
    dtR.Data{k}.ID = dtR.ID((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dtR.Data{k}.TX1 = dtR.X1((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dtR.Data{k}.TX2 = dtR.X2((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dtR.Data{k}.TY1 = dtR.Y1((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dtR.Data{k}.TY2 = dtR.Y2((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dtR.Data{k}.Vr = dtR.Vr((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dtR.Data{k}.Display = dtR.displayT((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    if strcmp(params.pType,'CircRDts')
        if mod(k,2) == 1
            % If natural gain
           dtR.Data{k}.VisRep = -dtR.VisualStimulusRep((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
        else
            % If reverse gain
            dtR.Data{k}.VisRep = dtR.VisualStimulusRep((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
        end
    end
    dtR.Data{k}.Vf = dtR.Vf((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dtR.Data{k}.Vs = dtR.Vs((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    [dtR.Data{k}.actState, dtR.Data{k}.ActBouts] = GetBouts(dtR.Data{k}.Vr, dtR.Data{k}.Vf, dtR.Data{k}.Vs);
    dtR.Data{k}.Frames = dtR.CameraFrames((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dtR.Data{k}.SoftTime = dtR.SoftwareTime((tTransitionTreadmill(k)+1):tTransitionTreadmill(k+1));
    dtR.Data{k}.FlickerSignalO = dtR.FlickerSignal(((tTransitionDAC(k)+1):tTransitionDAC(k+1)));
    dtR.Data{k}.MembranePotentialO = dtR.MembranePotential(((tTransitionDAC(k)+1):tTransitionDAC(k+1)));
    
    % Resize Flicker and Membrane potential to same time base as behavior
    tx = linspace(1,length(dtR.Data{k}.FlickerSignalO), length(dtR.Data{k}.FlickerSignalO));
    dtR.Data{k}.FlickerSignal = resample(dtR.Data{k}.FlickerSignalO,tx,...
        length(dtR.Data{k}.Vr)/length(dtR.Data{k}.FlickerSignalO), 'spline');
    dtR.Data{k}.MembranePotential = resample(dtR.Data{k}.MembranePotentialO,tx,...
        length(dtR.Data{k}.Vr)/length(dtR.Data{k}.MembranePotentialO), 'spline');

    % Resample Signals to the same time base 
    % Use intervals of 125ms for resampling
    delta = 500;
    ndiv = ceil((ddd(k)*params.SampRateTreadmill)/delta);
    rdiv = ceil(length(dtR.Data{k}.Vr)/ndiv);
    mdiv = ceil(length(dtR.Data{k}.MembranePotential)/ndiv);

    dtR.Data{k}.VrRes = [];
    dtR.Data{k}.DispRes = [];
    dtR.Data{k}.VisRepRes = [];
    dtR.Data{k}.VfRes = [];
    dtR.Data{k}.VsRes = [];
    dtR.Data{k}.acStRes = [];
    dtR.Data{k}.MPRes = [];
    dtR.Data{k}.FlickRes = [];
    if rdiv ~= 0
        for i = 1 : ndiv
            if (1+(i-1)*rdiv)<=length(dtR.Data{k}.Vr)
                if  i ~= ndiv && i*rdiv < length(dtR.Data{k}.Vr)
                    dtR.Data{k}.VrRes = vertcat(dtR.Data{k}.VrRes, resample(dtR.Data{k}.Vr((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    dtR.Data{k}.DispRes = vertcat(dtR.Data{k}.DispRes, resample(dtR.Data{k}.Display((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    dtR.Data{k}.VisRepRes = vertcat(dtR.Data{k}.VisRepRes, resample(dtR.Data{k}.VisRep((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    dtR.Data{k}.VfRes = vertcat(dtR.Data{k}.VfRes, resample(dtR.Data{k}.Vf((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    dtR.Data{k}.VsRes = vertcat(dtR.Data{k}.VsRes, resample(dtR.Data{k}.Vs((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    dtR.Data{k}.acStRes = vertcat(dtR.Data{k}.acStRes, resample(dtR.Data{k}.actState((1+(i-1)*rdiv):(i*rdiv)), ...
                        delta, rdiv, 0, 0));
                    dtR.Data{k}.MPRes = vertcat(dtR.Data{k}.MPRes, resample(dtR.Data{k}.MembranePotential((1+(i-1)*mdiv):(i*mdiv)), ...
                        delta, mdiv, 0, 0));
                    dtR.Data{k}.FlickRes = vertcat(dtR.Data{k}.FlickRes, resample(dtR.Data{k}.FlickerSignal((1+(i-1)*mdiv):(i*mdiv)), ...
                        delta, mdiv, 0, 0));
                else
                    auxrsp = resample(dtR.Data{k}.Vr((1+(i-1)*rdiv):end), delta, length(dtR.Data{k}.Vr((1+(i-1)*rdiv):end)), 0, 0);
                    dtR.Data{k}.VrRes = vertcat(dtR.Data{k}.VrRes, auxrsp(1:end));
                    auxrsp = resample(dtR.Data{k}.Display((1+(i-1)*rdiv):end), delta, length(dtR.Data{k}.Display((1+(i-1)*rdiv):end)), 0, 0);
                    dtR.Data{k}.DispRes = vertcat(dtR.Data{k}.DispRes, auxrsp(1:end));
                    auxrsp = resample(dtR.Data{k}.VisRep((1+(i-1)*rdiv):end), delta, length(dtR.Data{k}.VisRep((1+(i-1)*rdiv):end)), 0, 0);
                    dtR.Data{k}.VisRepRes = vertcat(dtR.Data{k}.VisRepRes, auxrsp(1:end));
                    auxrsp = resample(dtR.Data{k}.Vf((1+(i-1)*rdiv):end), delta, length(dtR.Data{k}.Vf((1+(i-1)*rdiv):end)), 0, 0);
                    dtR.Data{k}.VfRes = vertcat(dtR.Data{k}.VfRes, auxrsp(1:end));
                    auxrsp = resample(dtR.Data{k}.Vs((1+(i-1)*rdiv):end), delta, length(dtR.Data{k}.Vs((1+(i-1)*rdiv):end)), 0, 0);
                    dtR.Data{k}.VsRes = vertcat(dtR.Data{k}.VsRes, auxrsp(1:end));
                    auxrsp = resample(dtR.Data{k}.actState((1+(i-1)*rdiv):end), delta, length(dtR.Data{k}.actState((1+(i-1)*rdiv):end)), 0, 0);
                    dtR.Data{k}.acStRes = vertcat(dtR.Data{k}.acStRes, auxrsp(1:end));
                    auxrsp = resample(dtR.Data{k}.MembranePotential((1+(i-1)*mdiv):end), delta, length(dtR.Data{k}.MembranePotential((1+(i-1)*mdiv):end)), 0, 0);
                    dtR.Data{k}.MPRes = vertcat(dtR.Data{k}.MPRes,  auxrsp(1:end));
                    auxrsp = resample(dtR.Data{k}.FlickerSignal((1+(i-1)*mdiv):end), delta, length(dtR.Data{k}.FlickerSignal((1+(i-1)*mdiv):end)), 0, 0);
                    dtR.Data{k}.FlickRes = vertcat(dtR.Data{k}.FlickRes,  auxrsp(1:end));
                end
            end
        end
    end
end
% Save preprocess data
dtR.params = params;
dt = dtR;
disp(['Transitions Segmented: ' num2str(length(tTransitionDAC)-1)]);
save([path '\PreProcessedDataNoSR.mat'], 'dt','-v7.3')
disp('PreProcessedDataNoSR.mat file created')

end