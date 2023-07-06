function [params] = GetParamsEphys(protocol, PC, pathO)
% Parameters for pre processing Ephys data
if nargin == 1 % if protocol does not have subtrials
    PC = [];
end
if nargin == 2 % if we don't have to define output directory
    pathO = '';
end
params.pType = protocol; % store the used protocol

% Protocol types properties (if new protocols are to be added, add them as
% a new 'case' and define: ExperimentTime; ProtocolTransitionTime;
% Protocol; FlickerON; 
switch protocol
    case 'Replay'
        params.ExperimentTime = size(PC,2)*90; % Experiment duration
        params.ProtocolTransitionTime = [0 90:90:(size(PC,2)*90)]; % Transition times for subtrials
        params.Protocol = PC; % Sequence of subtrials
        params.FlickerON = 1; % Whether the photodiode to synchronize visual display and ephys is to be used
        params.pathO = pathO;
    case 'CircRDts'
        params.ExperimentTime = size(PC,2)*90; % Experiment duration
        params.ProtocolTransitionTime = [0 90:90:(size(PC,2)*90)]; % Transition times for subtrials
        params.Protocol = PC; % Sequence of subtrials
        params.FlickerON = 1; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'Flicker'
        params.ExperimentTime = 30; % Experiment duration
        params.ProtocolTransitionTime = [0 0.5:0.5:(60*0.5)]; % Transition times for subtrials
        params.Protocol = {'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', ...
            'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', ...
            'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', ...
            'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', ...
            'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', ...
            'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark', 'Light', 'Dark'};  % Sequence of subtrials
        params.FlickerON = 1; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'Dark1Min'
        params.ExperimentTime = 60; % Experiment duration
        params.ProtocolTransitionTime = [0 60]; % Transition times for subtrials
        params.Protocol = {'Dark'}; % Sequence of subtrials
        params.FlickerON = 0; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'Dark5Min'
        params.ExperimentTime = 300; % Experiment duration
        params.ProtocolTransitionTime = [0 300]; % Transition times for subtrials
        params.Protocol = {'Dark'}; % Sequence of subtrials
        params.FlickerON = 0; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'Dark30Min'
        params.ExperimentTime = 1800; % Experiment duration
        params.ProtocolTransitionTime = [0 1800]; % Transition times for subtrials
        params.Protocol = {'Dark'}; % Sequence of subtrials
        params.FlickerON = 0; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'MovShort'
        params.ExperimentTime = 53; % Experiment duration
        params.ProtocolTransitionTime = [0 2.5:2.5:(21*2.5)]; % Transition times for subtrials
        params.Protocol = {'Static', 'Right', 'Static', 'Left', 'Static', 'Right', 'Static', 'Left', ...
            'Static', 'Right', 'Static', 'Left', 'Static', 'Right', 'Static', 'Left', ...
            'Static', 'Right', 'Static', 'Left', 'Static'}; % Sequence of subtrials
        params.FlickerON = 1; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'RotDots'
        params.ExperimentTime = 256; % Experiment duration
        params.ProtocolTransitionTime = [0 4:4:(64*4)]; % Transition times for subtrials
        params.Protocol = {'1Static', '1Right', '1Static', '1Left', '2.5Static', '2.5Right', '2.5Static', '2.5Left', ...
            '5Static', '5Right', '5Static', '5Left', '10Static', '10Right', '10Static', '10Left', ...
            '1Static', '1Right', '1Static', '1Left', '2.5Static', '2.5Right', '2.5Static', '2.5Left', ...
            '5Static', '5Right', '5Static', '5Left', '10Static', '10Right', '10Static', '10Left', ...
            '1Static', '1Right', '1Static', '1Left', '2.5Static', '2.5Right', '2.5Static', '2.5Left', ...
            '5Static', '5Right', '5Static', '5Left', '10Static', '10Right', '10Static', '10Left', ...
            '1Static', '1Right', '1Static', '1Left', '2.5Static', '2.5Right', '2.5Static', '2.5Left', ...
            '5Static', '5Right', '5Static', '5Left', '10Static', '10Right', '10Static', '10Left'}; % Sequence of subtrials
        params.FlickerON = 1; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'MovMultiple'
        params.ExperimentTime = 256; % Experiment duration
        params.ProtocolTransitionTime = [0 4:4:(64*4)]; % Transition times for subtrials
        params.Protocol = {'16Static', '16Right', '16Static', '16Left', '32Static', '32Right', '32Static', '32Left', ...
            '64Static', '64Right', '64Static', '64Left', '128Static', '128Right', '128Static', '128Left', ...
            '16Static', '16Right', '16Static', '16Left', '32Static', '32Right', '32Static', '32Left', ...
            '64Static', '64Right', '64Static', '64Left', '128Static', '128Right', '128Static', '128Left', ...
            '16Static', '16Right', '16Static', '16Left', '32Static', '32Right', '32Static', '32Left', ...
            '64Static', '64Right', '64Static', '64Left', '128Static', '128Right', '128Static', '128Left', ...
            '16Static', '16Right', '16Static', '16Left', '32Static', '32Right', '32Static', '32Left', ...
            '64Static', '64Right', '64Static', '64Left', '128Static', '128Right', '128Static', '128Left'...
        }; % Sequence of subtrials
    params.FlickerON = 1; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'MovMultipleS'
        params.ExperimentTime = 160; % Experiment duration
        params.ProtocolTransitionTime = [0 2.5:2.5:(64*2.5)]; % Transition times for subtrials
        params.Protocol = {'16Static', '16Right', '16Static', '16Left', '32Static', '32Right', '32Static', '32Left', ...
            '64Static', '64Right', '64Static', '64Left', '128Static', '128Right', '128Static', '128Left', ...
            '16Static', '16Right', '16Static', '16Left', '32Static', '32Right', '32Static', '32Left', ...
            '64Static', '64Right', '64Static', '64Left', '128Static', '128Right', '128Static', '128Left', ...
            '16Static', '16Right', '16Static', '16Left', '32Static', '32Right', '32Static', '32Left', ...
            '64Static', '64Right', '64Static', '64Left', '128Static', '128Right', '128Static', '128Left', ...
            '16Static', '16Right', '16Static', '16Left', '32Static', '32Right', '32Static', '32Left', ...
            '64Static', '64Right', '64Static', '64Left', '128Static', '128Right', '128Static', '128Left'...
            }; % Sequence of subtrials
        params.FlickerON = 1; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'UniLeft'
        params.ExperimentTime = 53; % Experiment duration
        params.ProtocolTransitionTime = [0 2.5:2.5:(21*2.5)]; % Transition times for subtrials
        params.Protocol = {'Static', 'BTF', 'Static', 'FTB', 'Static', 'BTF', 'Static', 'FTB', ...
            'Static', 'BTF', 'Static', 'FTB', 'Static', 'BTF', 'Static', 'FTB', ...
            'Static', 'BTF', 'Static', 'FTB', 'Static'}; % Sequence of subtrials
        params.FlickerON = 1; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'UniRight'
        params.ExperimentTime = 53; % Experiment duration
        params.ProtocolTransitionTime = [0 2.5:2.5:(21*2.5)]; % Transition times for subtrials
        params.Protocol = {'Static', 'FTB', 'Static', 'BTF', 'Static', 'FTB', 'Static', 'BTF', ...
            'Static', 'FTB', 'Static', 'BTF', 'Static', 'FTB', 'Static', 'BTF', ...
            'Static', 'FTB', 'Static', 'BTF', 'Static'}; % Sequence of subtrials
        params.FlickerON = 1; % Whether the photodiode to synchronize visual display and ephys is to be used
    case 'MovMultiple3'
        params.ExperimentTime = 70; % Experiment duration
        params.ProtocolTransitionTime = [0 10:10:70]; % Transition times for subtrials
        params.Protocol = {'Stat', 'Mov', 'Stat', 'Mov', 'Stat', 'Mov', 'Stat'}; % Sequence of subtrials
    params.FlickerON = 1; % Whether the photodiode to synchronize visual display and ephys is to be used
    otherwise
        params.ExperimentTime = 0;
        params.ProtocolTransitionTime = 0;
        params.Protocol = {};
end
% General sampling parameters
params.SampRateDAC = 10000;
params.SampRateTreadmill = 4000;
params.PhotodiodeFlickRate = 400;
params.TargetFrameRate = 500; % Hz
params.SmoothParameterBehavior = 2000;

% General calibration parameters
params.CalibRot = 0.449; % º/tick   positive -> Leftward?
params.CalibForw = -0.0257; % mm/tick
end

