clear
clc
% Add directories of closed loop trials to pre-process
path{001} = '\';

% Loop for every closed loop trial
for i = 1 : length(path)
    protocol = ''; % Protocol Type (see Params)
    PC = {}; % Protocol subtrial sequence
    if ~isempty(path{i})
        [params] = GetParamsEphys(protocol, PC, path{i}); % Load parameters
        [dt] = LoadData(path{i}); % Load all time series data
        [params, dt] = ProcessClockAndTriggers(params, dt, params.FlickerON); % Process triggers, measure lags and align different aquisitions
        [dt] = ProcessRawDataNoSR(params, dt, path{i}); % Crop the raw traces into subtrials and pre-process ball data into behavior traces
    end
    disp([num2str(i) ' Done.'])
end


%%
clear
clc
% Add directories of replay and respective closed loop trials to pre-process
pathR{001} ='\';
pathV{001} ='\';

% Loop for every replay trial
for i = 1 : length(pathR)
    protocol = ''; % Protocol Type (see Params)
    PC = {}; % Protocol subtrial sequence
    if ~isempty(pathR{i})
        [params] = GetParamsEphys(protocol, PC, pathR{i});  % Load parameters
        [dt] = LoadData(pathR{i}); % Load all time series data from replay trial
        [params, dt] = ProcessClockAndTriggers(params, dt, params.FlickerON); % Process triggers, measure lags and align different aquisitions
        dtCL = load([pathV{i} 'PreProcessedDataNoSR.mat']); % Load pre process respective closed loop trial
        dtCL = dtCL.dt;
        [dt] = ProcessRawDataNoSRReplay(params, dt, dtCL, pathR{i});  % Crop the raw traces into subtrials, pre-process ball data into behavior traces and add the aligned replayed visual stimulus
    end
end




















