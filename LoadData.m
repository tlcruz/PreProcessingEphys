function [dt] = LoadData(path)
% Load the treadmill time series
f1 = load([path 'f1.txt']);
% Load the ephys series
Data = fileread([path 'f2.txt']);
Data = strrep(Data, ',', '.'); % Replace , by .
FID = fopen([path 'f2.txt'], 'w');
fwrite(FID, Data, 'char');
fclose(FID);
clearvars Data FID
f2 = load([path 'f2.txt']); % Reload

% Attribute different time series to labeled properties
dt.ID = f1(:,1); % frame ID
dt.TClock = f1(:,1); % treadmill clock signal
dt.X1 = f1(:,3); % horizontal channel treadmill camera 1
dt.X2 = f1(:,5); % horizontal channel treadmill camera 2
dt.Y1 = f1(:,4); % vertical channel treadmill camera 1
dt.Y2 = f1(:,6); % vertical channel treadmill camera 2
dt.SoftwareTime = f1(:,7); % software clock 
dt.CameraFrames = f1(:,8); % frame number recording camera
dt.MembranePotential = 100 * f2(:,1); % membrane potential
dt.FlickerSignal = f2(:,2); % flicker signal from the photodiode
dt.TreadmillTrigger = f2(:,3); % treadmill acquisition trigger

disp([path '...  Data Loaded Successfully'])
end