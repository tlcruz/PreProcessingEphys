function GuiPlotRawDataCLCurateSR(path)
%% Get Directory and Individual Fly Folders
if nargin ==0
    path = uigetdir('C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\DATA\TOMÁS\Ephys VR\Data\H2 Recordings');
    path = [path '\'];
end
flies = dir([path, '\*']);
flies = flies(3:end);
flies(~[flies.isdir]) = [];
currentFly = 1;
currentTrial = 1;
currentPart = 1;

%% Load DataFile Current Fly
dt = load([path flies(currentFly).name '\Trial 0' num2str(currentTrial) '\PreProcessedDataNoSR.mat']);
dt = dt.dt;
%% Figure Properties
f = figure('Visible', 'off', 'Position', [360, 300, 500, 320]);
f.Name = 'NavigateVideo';
movegui(f, 'center')
f.Visible = 'on';
ha = axes('Units', 'pixels', 'Position', [50, 50, 430, 200]);
ha.Units = 'normalized';
axes(ha);
axis off;

%% Define GUI elements
% Fly Number
NFlyplus = uicontrol('Style', 'pushbutton', 'String', '+',...
    'Position', [40, 291, 9, 14], 'Callback', {@NFlyplus_Callback});
NFlyplus.Units = 'normalized';
NFlyminus = uicontrol('Style', 'pushbutton', 'String', '-',...
    'Position', [40, 275, 9, 14], 'Callback', {@NFlyminus_Callback});
NFlyminus.Units = 'normalized';
NFlyT = uicontrol(f,'Style','edit','String','FlyNo',...
    'Position',[20, 275, 20, 30], 'Callback', {@NFlyT_CallBack});
NFlyT.Units = 'normalized';

% Trial Number
NTrialplus = uicontrol('Style', 'pushbutton', 'String', '+',...
    'Position', [70, 298, 9, 7], 'Callback', {@NTrialplus_Callback});
NTrialplus.Units = 'normalized';
NTrialminus = uicontrol('Style', 'pushbutton', 'String', '-',...
    'Position', [70, 291, 9, 7], 'Callback', {@NTrialminus_Callback});
NTrialminus.Units = 'normalized';
NTrialT = uicontrol(f,'Style','edit','String','TrialNo',...
    'Position',[50, 291, 20, 14], 'Callback', {@NTrialT_CallBack});
NTrialT.Units = 'normalized';

% Trial part
NPartplus = uicontrol('Style', 'pushbutton', 'String', '+',...
    'Position', [70, 282, 9, 7],'Callback', {@NPartplus_Callback});
NPartplus.Units = 'normalized';
NPartminus = uicontrol('Style', 'pushbutton', 'String', '-', ...
    'Position', [70, 275, 9, 7],'Callback', {@NPartminus_Callback});
NPartminus.Units = 'normalized';
NPartT = uicontrol(f,'Style','edit','String','PartNo',...
    'Position',[50, 275, 20, 14],'Callback', {@NPartT_CallBack});
NPartT.Units = 'normalized';

% Fly number
FlyNr = uicontrol(f,'Style','text','String','Fly Nr.',...
    'Position',[20, 305, 50, 10]);
FlyNr.Units = 'normalized';
PlotFly = uicontrol('Style', 'pushbutton', 'String', 'Save',...
    'Position', [350, 280, 20, 10], 'Callback', {@PlotFly_Callback});
PlotFly.Units = 'normalized';

% Voltage multiplier
AVMPText = uicontrol(f,'Style','text','String','MPMultiplier',...
    'Position',[100, 295, 20, 10]);
AVMPText.Units = 'normalized';

AVMP = uicontrol(f,'Style','edit','String','Multiplier',...
    'Position',[100, 280, 20, 10],'Callback', {@Multiplier_CallBack});
AVMP.Units = 'normalized';

% Filtered threshold
SRText = uicontrol(f,'Style', 'pushbutton', 'String','SetSpkThr',...
    'Position',[130, 295, 20, 10], 'Callback', {@SetThreshold_Callback});
SRText.Units = 'normalized';

SRMP = uicontrol(f,'Style','edit','String','SpkThr',...
    'Position',[130, 280, 20, 10], 'Callback', {@SpikeThreshold_CallBack});
SRMP.Units = 'normalized';

% Use high pass filter
CCText = uicontrol(f,'Style', 'pushbutton', 'String','UseHPFilter',...
    'Position',[200, 295, 20, 10], 'Callback', {@UseHPFilter_Callback});
CCText.Units = 'normalized';

% High pass filter cutoff
CCT = uicontrol(f,'Style','edit','String','HP Filter Freq.',...
    'Position',[200, 280, 20, 10], 'Callback', {@HPFilterTreshold_Callback});
CCT.Units = 'normalized';

% Gaussian window window
MuText = uicontrol(f,'Style','text','String','Gaussian Window',...
    'Position',[230, 295, 20, 10]);
MuText.Units = 'normalized';

muT = uicontrol(f,'Style','edit','String','Window',...
    'Position',[230, 280, 20, 10],'Callback', {@Window_CallBack});
muT.Units = 'normalized';

% Use 50Hz notch filter
SigText = uicontrol(f,'Style','pushbutton','String','Use50HzFilter',...
    'Position',[260, 295, 20, 10], 'Callback', {@Use50HzFilter_Callback});
SigText.Units = 'normalized';

% Check if data is not usable
LText = uicontrol(f,'Style','pushbutton','String','NonUsableData',...
    'Position',[300, 295, 20, 10], 'Callback', {@IsUsable_Callback});
LText.Units = 'normalized';

% Display subtrials
ww = ceil(430/length(dt.params.Protocol));
for n = 1 : length(dt.params.Protocol)
    if n == floor(currentPart)
        nseq{n} = uicontrol('Style','text','position',[ww*n, 290-27, ww, 10],'String', dt.params.Protocol{n}, 'BackgroundColor', 'r');
    else
        nseq{n} = uicontrol('Style','text','position',[ww*n, 290-27, ww, 10],'String', dt.params.Protocol{n}, 'BackgroundColor', 'w');
    end
    nseq{n}.Units = 'normalized';
end

% Initialize transformations to calculate spike rate
pTransform.MPMultiplier = 1*ones(length(dt.params.Protocol),1);
pTransform.noiseWallFilter = 1*ones(length(dt.params.Protocol),1);
pTransform.HPFilter = 1*ones(length(dt.params.Protocol),1);
pTransform.HPcutoff = 350*ones(length(dt.params.Protocol),1);
pTransform.SetSpkThr = 0*ones(length(dt.params.Protocol),1);
pTransform.SpkThr = 2.5*ones(length(dt.params.Protocol),1);
pTransform.GaussWindow = 800*ones(length(dt.params.Protocol),1);
pTransform.UsableData = 1*ones(length(dt.params.Protocol),1);

%% CallBacks
% Store the spike rate calculated using the defined transformations
    function PlotFly_Callback(source, eventdata)
        if currentFly <= length(flies)
            disp('Generating PreProcessedDataWithSR.mat file')
            % Apply multiplier
            MPOri = dt.MembranePotential; 
            if pTransform.MPMultiplier(1) ~= 1
                MPOri = MPOri/pTransform.MPMultiplier(1);
            end
            % Apply notch filter 50Hz
            if pTransform.noiseWallFilter(1) == 1 
                d3 = designfilt('bandstopiir','FilterOrder',10, ...
                    'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
                    'SampleRate',dt.params.SampRateDAC);
                MPOri2 = filtfilt(d3, MPOri);
            end
            % Apply high pass filter
            if pTransform.HPFilter(1) == 1 
                d1 = designfilt('highpassiir','FilterOrder',12, ...
                    'HalfPowerFrequency',pTransform.HPcutoff(1)/dt.params.SampRateDAC,'DesignMethod','butter');
                if pTransform.noiseWallFilter(1) ==1
                    MPOri2 = filtfilt(d1, MPOri2);
                else
                    MPOri2 = filtfilt(d1, MPOri);
                end
            end
            % find peaks in the filtered trace
            [~, locs, w, prom]= findpeaks(MPOri2);
            % Keep peaks shorter than 3ms
            locs = locs(w<30);
            prom = prom(w<30);
            % Get the histogram or peak prominences and find the threshold 
            % between noise and spikes
            promCnt = min(prom) : 0.1: max(prom);
            aux = smooth(hist(prom, promCnt));
            [~,lc, ~, p] = findpeaks(-aux);
            [~,IR] = max(p);
            spkThr = promCnt(lc(max(IR)));
            % Manually set threshold
            if pTransform.SetSpkThr(1) == 1
                spkThr = pTransform.SpkThr(1);
            end
            if isempty(spkThr)
                spkThr = 1;
            end
            % Keep the spikes over the defined threshold
            locs = locs(prom>spkThr);
            spikes = zeros(length(MPOri),1);
            if ~isempty(spkThr)
                spikes(locs) = 1;
            end
            % Convolve with a gaussian window
            w = gausswin(pTransform.GaussWindow(1));
            dt.SROut = conv(spikes,w,'same');
            % Do spike rate detection per subtrial
            for k = 1 : length(dt.Data)
                MPOri = dt.Data{k}.MembranePotentialO;
                VrOut = dt.Data{k}.Vr;
                % Apply multiplier
                if pTransform.MPMultiplier(k) ~= 1
                    MPOri = MPOri/pTransform.MPMultiplier(k);
                end
                % Apply notch filter 50Hz
                if pTransform.noiseWallFilter(k) == 1
                    d3 = designfilt('bandstopiir','FilterOrder',10, ...
                        'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
                        'SampleRate',dt.params.SampRateDAC);
                    MPOri2 = filtfilt(d3, MPOri);
                end
                % Apply high pass filter
                if pTransform.HPFilter(k) == 1
                    d1 = designfilt('highpassiir','FilterOrder',12, ...
                        'HalfPowerFrequency',pTransform.HPcutoff(k)/dt.params.SampRateDAC,'DesignMethod','butter');
                    if pTransform.noiseWallFilter(k) ==1
                        MPOri2 = filtfilt(d1, MPOri2);
                    else
                        MPOri2 = filtfilt(d1, MPOri);
                    end
                end
                % find peaks in the filtered trace
                [~, locs, w, prom]= findpeaks(MPOri2);
                % Keep peaks shorter than 3ms
                locs = locs(w<30);
                prom = prom(w<30);
                % Get the histogram or peak prominences and find the threshold 
                % between noise and spikes
                promCnt = min(prom) : 0.1: max(prom);
                aux = smooth(hist(prom, promCnt));
                [~,lc, ~, p] = findpeaks(-aux);
                [~,IR] = max(p);
                spkThr = promCnt(lc(max(IR)));
                % Set the threshold manually
                if pTransform.SetSpkThr(k) == 1
                    spkThr = pTransform.SpkThr(k);
                end
                if isempty(spkThr)
                    spkThr = 1;
                end
                pTransform.SetSpkThr(k)
                % Keep the spikes over the defined threshold
                locs = locs(prom>spkThr);
                spikes = zeros(length(MPOri),1);
                if ~isempty(spkThr)
                    spikes(locs) = 1;
                end
                dt.Data{k}.spikesO = spikes;
                % Convolve with a gaussian window
                w = dt.params.SampRateDAC/pTransform.GaussWindow(k) * gausswin(pTransform.GaussWindow(k));
                SROut = conv(spikes,w,'same');
                % Resample calculated spike rate into the behavior time base
                dt.Data{k}.SRO = SROut;
                tx = linspace(1,length(SROut), length(SROut));
                auxSpikes = resample(spikes,tx,length(dt.Data{k}.Vr)/length(spikes), 'pchip');
                auxSpikes(auxSpikes<0.2) = 0;
                auxSpikes(auxSpikes>0) = 1;
                auxSpikes = vertcat(0,diff(auxSpikes));
                auxSpikes(auxSpikes>-0.8) = 0;
                auxSpikes(auxSpikes<=-0.8) = 1;
                dt.Data{k}.spikes = auxSpikes;
                tx = linspace(1,length(SROut), length(SROut));
                dt.Data{k}.SR = resample(SROut,tx,length(VrOut)/length(SROut), 'spline');
            end
            % Save new file including the curated spike rate
            dt.SrPars = pTransform;
            save([path flies(currentFly).name '\Trial 0' num2str(currentTrial) '\PreProcessedDataWithSR.mat'], 'dt','-v7.3')
            disp('PreProcessedDataWithSR.mat file created')
        end
    end
    % Function that handles the multiplier input
    function Multiplier_CallBack(hObject, eventdata, handles)
         input = str2double(get(hObject,'String'));
         if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
         else
             pTransform.MPMultiplier(currentPart) = input;
             updateAxes(); % refresh the GUI display
         end
    end
    % Function that handles the yes/no threshold input
    function SetThreshold_Callback(source, eventdata)
        if pTransform.SetSpkThr(currentPart) == 1
            pTransform.SetSpkThr(currentPart) = 0;
        else
            pTransform.SetSpkThr(currentPart) = 1;
        end
        updateAxes(); % refresh the GUI display
    end
    % Function that handles the threshold input
    function SpikeThreshold_CallBack(hObject, eventdata, handles)
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            pTransform.SpkThr(currentPart) = input;
            if pTransform.SetSpkThr(currentPart) == 1
                pTransform.SetSpkThr(currentPart) = 0;
            else
                pTransform.SetSpkThr(currentPart) = 1;
            end
            updateAxes(); % refresh the GUI display
        end
    end
    % Function that handles the yes/no high pass filter input
    function UseHPFilter_Callback(source, eventdata)
        if pTransform.HPFilter(currentPart) == 1
            pTransform.HPFilter(currentPart) = 0;
        else
            pTransform.HPFilter(currentPart) = 1;
        end
        updateAxes(); % refresh the GUI display
    end
    % Function that handles the high pass filter input
    function HPFilterTreshold_Callback(hObject, eventdata, handles)
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            pTransform.HPcutoff(currentPart) = input;
            updateAxes(); % refresh the GUI display
        end
    end
    % Function that handles the gaussian window input
    function Window_CallBack(hObject, eventdata, handles)
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            pTransform.GaussWindow(currentPart) = input;
            updateAxes(); % refresh the GUI display
        end
    end
    % Function that handles the 50 Hz notch filter input
    function Use50HzFilter_Callback(source, eventdata)
        if pTransform.noiseWallFilter(currentPart) == 1
            pTransform.noiseWallFilter(currentPart) = 0;
        else
            pTransform.noiseWallFilter(currentPart) = 1;
        end
        updateAxes(); % refresh the GUI display
    end
    % Function that handles the data usability input
    function IsUsable_Callback(source, eventdata)
        if pTransform.UsableData(currentPart) == 1
            pTransform.UsableData(currentPart) = 0;
        else
            pTransform.UsableData(currentPart) = 1;
        end
        updateAxes(); % refresh the GUI display
    end
    % Function that handles the fly ID increase input
    function NFlyplus_Callback(source, eventdata)
        if currentFly <= length(flies)
            currentFly = currentFly + 1;
            dt = load([path flies(currentFly).name '\Trial 0' num2str(currentTrial) '\PreProcessedDataNoSR.mat']);
            dt = dt.dt;
            updateAxes();
        end
    end
    % Function that handles the fly ID decrease input
    function NFlyminus_Callback(source, eventdata)
        if currentFly > 1
            currentFly = currentFly - 1;
            dt = load([path flies(currentFly).name '\Trial 0' num2str(currentTrial) '\PreProcessedDataNoSR.mat']);
            dt = dt.dt;
            updateAxes(); % refresh the GUI display
        end
    end
    % Function that handles the fly ID manual input
    function NFlyT_CallBack(hObject, eventdata, handles)
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            input = floor(input);
            if input <= length(flies) && input >= 1
                currentFly = input;
                dt = load([path flies(currentFly).name '\Trial 0' num2str(currentTrial) '\PreProcessedDataNoSR.mat']);
                dt = dt.dt;
                updateAxes(); % refresh the GUI display
            else
                currentFly = 1;
                dt = load([path flies(currentFly).name '\Trial 0' num2str(currentTrial) '\PreProcessedDataNoSR.mat']);
                dt = dt.dt;
                updateAxes(); % refresh the GUI display
            end
        end
    end
    % Function that handles the trial ID increase input
    function NTrialplus_Callback(source, eventdata)
        if currentTrial <= length(flies)
            currentTrial = currentTrial + 1;
            dt = load([path flies(currentFly).name '\Trial 0' num2str(currentTrial) '\PreProcessedDataNoSR.mat']);
            dt = dt.dt;
            updateAxes(); % refresh the GUI display
        end
    end
    % Function that handles the trial ID decrease input
    function NTrialminus_Callback(source, eventdata)
        if currentTrial > 1
            currentTrial = currentTrial - 1;
            dt = load([path flies(currentFly).name '\Trial 0' num2str(currentTrial) '\PreProcessedDataNoSR.mat']);
            dt = dt.dt;
            updateAxes(); % refresh the GUI display
        end
    end
    % Function that handles the trial ID manual input
    function NTrialT_CallBack(hObject, eventdata, handles)
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            input = floor(input);
            if input <= length(flies) && input >= 1
                currentTrial = input;
                updateAxes(); % refresh the GUI display
            else
                currentTrial = 1;
                updateAxes(); % refresh the GUI display
            end
        end
    end
    % Function that handles the subtrial ID increase input
    function NPartplus_Callback(source, eventdata)
        if currentPart <= length(dt.params.Protocol)-2
            currentPart = currentPart + 1;
            updateAxes(); % refresh the GUI display
        end
    end
    % Function that handles the subtrial ID decrease input
    function NPartminus_Callback(source, eventdata)
        if currentPart > 1
            currentPart = currentPart - 1;
            updateAxes(); % refresh the GUI display
        end
    end
    % Function that handles the subtrial ID manual input
    function NPartT_CallBack(hObject, eventdata, handles)
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            input = floor(input);
            if input <= length(dt.params.Protocol) && input >= 1
                currentPart = input;
                updateAxes(); % refresh the GUI display
            else
                currentPart = 1;
                updateAxes(); % refresh the GUI display
            end
        end
    end

%% UpdateAxes
    function updateAxes()
        % Update labels to represent the data being displayed
        FlyNr.String = flies(currentFly).name;
        NPartT.String = num2str(currentPart);
        NTrialT.String = num2str(currentTrial);
        NFlyT.String = num2str(currentFly);
        for k = 1 : length(dt.params.Protocol)
            if k == floor(currentPart)
                set(nseq{k},'BackgroundColor',[1 0 0]);
            else
                set(nseq{k},'BackgroundColor',[1 1 1]);
            end
        end
        
        % Load the time series of the subtrial to display
        MPOri = dt.Data{currentPart}.MembranePotentialO;
        VrOut = dt.Data{currentPart}.Vr;
        VfOut = dt.Data{currentPart}.Vf;
        VsOut = dt.Data{currentPart}.Vs;
        VisOut = dt.Data{currentPart}.VisCL;
        % Apply voltage multiplier
        if pTransform.MPMultiplier(currentPart) ~= 1
            MPOri = MPOri/pTransform.MPMultiplier(currentPart);
        end
        % Apply 50Hz notch filter
        if pTransform.noiseWallFilter(currentPart) == 1
            d3 = designfilt('bandstopiir','FilterOrder',10, ...
                'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
                'SampleRate',dt.params.SampRateDAC);
            MPOri2 = filtfilt(d3, MPOri);
        end
        % Apply high pass filter
        if pTransform.HPFilter(currentPart) == 1
            d1 = designfilt('highpassiir','FilterOrder',12, ...
                'HalfPowerFrequency',pTransform.HPcutoff(currentPart)/dt.params.SampRateDAC,'DesignMethod','butter');
            if pTransform.noiseWallFilter(currentPart) ==1
                MPOri2 = filtfilt(d1, MPOri2);
            else
                MPOri2 = filtfilt(d1, MPOri);
            end
        end
        Tt = (1:length(MPOri))/dt.params.SampRateDAC;
        % find peaks in the filtered trace
        [~, locs, w, prom]= findpeaks(MPOri2);
        % Keep peaks shorter than 3ms
        locs = locs(w<30);
        prom = prom(w<30);
        % Get the histogram or peak prominences and find the threshold
        promCnt = min(prom) : 0.1: max(prom);
        aux = smooth(hist(prom, promCnt));
        if length(aux) > 3
            [~,lc, ~, p] = findpeaks(-aux);
            [~,IR] = max(p);
            spkThr = promCnt(lc(max(IR)));
        else
            spkThr = 1;
        end
        % Set the threshold manually
        if pTransform.SetSpkThr(currentPart) == 1
            spkThr = pTransform.SpkThr(currentPart);
        end
        if isempty(spkThr)
            spkThr = 1;
        end
        % Keep the spikes over the defined threshold
        locs = locs(prom>spkThr);
        spikes = zeros(length(MPOri),1);
        if ~isempty(spkThr)
            spikes(locs) = 1;
        end
        % Convolve with a gaussian window
        w = dt.params.SampRateDAC/pTransform.GaussWindow(currentPart) * gausswin(pTransform.GaussWindow(currentPart));
        SROut = conv(spikes,w,'same');
        % Resample calculated spike rate into the behavior time base
        tx = linspace(1,length(SROut), length(SROut));
        auxSpikes = resample(spikes,tx,length(VrOut)/length(spikes), 'pchip');
        auxSpikes(auxSpikes<0.2) = 0;
        auxSpikes(auxSpikes>0) = 1;
        auxSpikes = vertcat(0,diff(auxSpikes));
        auxSpikes(auxSpikes>-0.8) = 0;
        auxSpikes(auxSpikes<=-0.8) = 1;
        SRRes = resample(SROut,tx,length(VrOut)/length(SROut), 'spline');
        
        % Plot the time series
        Ttread = (1:length(SRRes))/dt.params.SampRateTreadmill;
        subplot(7,4,[5 6 7 9 10 11 13 14 15 17 18 19 21 22 23 25 26 27])
        cla reset;
        hold on
        plot(Tt, MPOri2+(nanmean(MPOri)-nanmean(MPOri2)), 'b', 'linewidth', 0.1) % Plot filtered membrane potential
        plot(Tt, -20+0.1*SROut, 'k', 'linewidth', 0.5)  % Plot spike rate
        plot(Ttread, -20+auxSpikes, 'r', 'linewidth', 0.5)  % Plot spike times
        scatter(Tt(locs), MPOri(locs), 50, 'r') % Plot spikes
        plot(Tt, MPOri, 'k', 'linewidth', 0.2) % Plot membrane potential
        plot(Ttread, -50 + 0.1*VrOut, 'Color', [0,0,1], 'linewidth', 0.5) % Plot angular velocity
        plot(Ttread, -50 + 0.1*VisOut, 'Color', [1,0,0], 'linewidth', 0.5) % Plot forward velocity
        axis([1 360000/dt.params.SampRateTreadmill -65 -10])
        
        % Plot the membrane potential spectrum
        subplot(7,4,[8 12])
        cla reset;
        hold on
        Fs = dt.params.SampRateDAC;
        % Calculate the FFT
        L = length(MPOri);
        Y = fft(MPOri);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        plot(f,log(P1), 'k') % Plot the spectrum
        axis([0 500 -15 5])
        ylabel('log|FFT(f)|')
        
        % Plot the histogram of spike prominences
        subplot(7,4,[16 20])
        cla reset;
        hold on
        plot(promCnt, log(aux), 'k')
        plot([spkThr spkThr], [-5 20], 'r')
        ylabel('Spike prominence histogram')
        
        % Plot the covariance between spike rate and angular velocity
        minL = 2*4000;
        % Isolate moments when the fly or the visual environment are active
        [~, ActBoutsR] = GetBouts(VrOut, VfOut, VsOut, 1000, 500, 5, 0.1);
        [~, ActBoutsV] = GetBouts(VisOut, VfOut, VsOut, 1000, 500, 5, 0.1);
        cvsR = [];
        for ii = 1 : length(ActBoutsR)
            if length(ActBoutsR{ii}) > 2000
                % Calculate cross covariance during  fly movement
                [cM,lagsM] = xcov(VrOut(ActBoutsR{ii}), SRRes(ActBoutsR{ii}), minL/2, 'coeff');
                cvsR = horzcat(cvsR, cM);
            end
        end
        
        cvsV = [];
        minL = 2*4000;
        for ii = 1 : length(ActBoutsV)
            if length(ActBoutsV{ii}) > 2000
                % Calculate cross covariance during visual environment movement
                [cV,lagsV] = xcov(VisOut(ActBoutsV{ii}), SRRes(ActBoutsV{ii}), minL/2, 'coeff');
                cvsV = horzcat(cvsV, cV);
            end
        end
        
        %  Plot the covariances and the timing of the covariance peak
        subplot(7,4,[24 28])
        cla reset;
        hold on
        if ~isempty(cvsR)
            [~,IR] = max(abs(mean(cvsR,2)));
            [~,IV] = max(abs(mean(cvsV,2)));
            plot([-1 1], [0 0], '--k', 'linewidth', 0.5)
            plot([0 0],[-1 1], '--k', 'linewidth', 0.5)
            plot([lagsM(IR) lagsM(IR)]/dt.params.SampRateTreadmill, [-1 1], '--b', 'linewidth', 0.5)
            plot(lagsM/dt.params.SampRateTreadmill, mean(cvsR,2), 'b', 'linewidth', 3)
            plot(lagsM/dt.params.SampRateTreadmill, mean(cvsR,2)+std(cvsR,1,2), 'b', 'linewidth', 1)
            plot(lagsM/dt.params.SampRateTreadmill, mean(cvsR,2)-std(cvsR,1,2), 'b', 'linewidth', 1)
            plot([lagsV(IV) lagsV(IV)]/dt.params.SampRateTreadmill, [-1 1], '--r', 'linewidth', 0.5)
            plot(lagsV/dt.params.SampRateTreadmill, mean(cvsV,2), 'r', 'linewidth', 3)
            plot(lagsV/dt.params.SampRateTreadmill, mean(cvsV,2)+std(cvsV,1,2), 'r', 'linewidth', 1)
            plot(lagsV/dt.params.SampRateTreadmill, mean(cvsV,2)-std(cvsV,1,2), 'r', 'linewidth', 1)
            axis([-minL/(2*dt.params.SampRateTreadmill) minL/(2*dt.params.SampRateTreadmill) -1 1])
        end
    end
end
