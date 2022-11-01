%% activPAL tri-monitor joint angle determination for sedentary bouts
% Requires events and full data file for each monitor. Written with r2022a
% Author: W. Seth Daley, Oct 7th 2022
% Program loads event and full data files as exported by ActivPAL 3 or 4 

% Copyright 2022 W Seth Daley
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%    http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License. 



%% Import data files
% Two data files are required for each activPAL, the events file, and the
% file with the full, uncompressed acceleration data. The user will be
% prompted to select the 6 files for the associated participant. Without
% some kind of file naming convention this process cannot be automated. 

clc
clear all

% Get user to select the Thigh files
[org.thigh.EventsFileName,org.thigh.FileLoc] = uigetfile({'*.csv',...
    'CSV Files (*.csv)';'*.*','All Files'},...
    'Select the Thigh Events file');

[org.thigh.FullFileName,~] = uigetfile({'*.csv',...
    'CSV Files (*.csv)';'*.*','All Files'},...
    'Select the Thigh full data file',org.thigh.FileLoc);

% Get user to select the Tibia files
[org.tibia.EventsFileName,org.tibia.FileLoc] = uigetfile({'*.csv',...
    'CSV Files (*.csv)';'*.*','All Files'},...
    'Select the Tibia Events file',org.thigh.FileLoc);

[org.tibia.FullFileName,~] = uigetfile({'*.csv',...
    'CSV Files (*.csv)';'*.*','All Files'},...
    'Select the Tibia full data file',org.tibia.FileLoc);

% Get user to select the Torso events file
[org.torso.EventsFileName,org.torso.FileLoc] = uigetfile({'*.csv',...
    'CSV Files (*.csv)';'*.*','All Files'},...
    'Select the Torso Events file',org.tibia.FileLoc);

[org.torso.FullFileName,~] = uigetfile({'*.csv',...
    'CSV Files (*.csv)';'*.*','All Files'},...
    'Select the Torso full data file',org.torso.FileLoc);

% Set a variable with the names of the three activpal placements
org.segmentNames = {'torso','thigh','tibia'};

% Create a wait bar for the file loading
f = waitbar(0, 'Getting ready to load files');
n = 6;
c = 0;

% Loop over the number of segments
for i=1:3
    
    % load the segment name
    name = org.segmentNames{i};
    
    % Advance the waitbar counter
    c = c+1;

    % Create a string with the location and name of the events file
    temp.(name).EventName = strcat(org.(name).FileLoc,org.(name).EventsFileName);

    % Load the events data
    raw.eventData.(name) = readmatrix(temp.(name).EventName);

    % Display the waitbar
    waitbar(c/n, f, sprintf('Loading files: %d %%', floor(c/n*100)));

    % Depending on the output, there are sometimes initial rows of zeros
    % that should be removed
    while raw.eventData.(name)(1,1) == 0
        raw.eventData.(name)(1,:) = [];
    end

    % Advance the waitbar
    c = c+1;

    % Create a string with the location and name of the full acceleration file
    temp.(name).FullDataName = strcat(org.(name).FileLoc,org.(name).FullFileName);

    % Load the full acceleration file
    raw.data.(name) = readmatrix(temp.(name).FullDataName);

    % Depending on the version of software used to download the data from
    % the activPAL there may be an additional column of unneeded data in
    % the acceleration file that should be removed
    % raw.data.(name)(:,2) = [];

    % Display the waitbar
    waitbar(c/n, f, sprintf('Loading files: %d %%', floor(c/n*100)));

    % Remove variable used within the loop
    clear name
end

close(f)
clear p n f c

% clear temporary variables
clear temp* i

% set variables (sample frequency, AP version, and the frequency used to
% filter acceleration data)
x = inputdlg({'What was the sample frequency?',...
    'Which version of activPAL was used?','Desired filter frequency','Downsample cutoff (h)'},...
    'Input',[1 20],{'20','3','0.18','2'});
org.sampleFreq = str2double(x(1));
org.APversion = str2double(x(2));
org.filterFreq = str2double(x(3));
org.downsampleCutoff = str2double(x(4));
clear x

% Check that the AP version is supported
if (org.APversion ~= 3) && (org.APversion ~= 4)
    error('That version of Activpal is not supported')
end

%% Find sedentary bouts

% Using thigh data to determine sedentary bouts as per prior research 

% Initialize a variable to store the time codes for sedentary events
org.eventTimeCodes = [];
org.eventDuration = [];

% Loop over the data with position IDs
for i=1:length(raw.eventData.thigh)-1
    
    % Find points where the position code is 0 (sedentary)
    if (raw.eventData.thigh(i,4) == 0)
        
        % Store the time associated with that index in a temporary variable
        tempTime = raw.eventData.thigh(i,1);
        tempDuration = raw.eventData.thigh(i,3);
        
        % Append the temporary variable to the list of event times
        org.eventTimeCodes = [org.eventTimeCodes tempTime];
        org.eventDuration = [org.eventDuration tempDuration];

    end
    clear temp*
end
clear i

% Record the number of events
org.numEvents = length(org.eventTimeCodes);

%% Crop sedentary events from data

% Loop over the segments
for i=1:3
   
    % load the segment name as a temporary variable
    name = org.segmentNames{i};

    % Loop over the number of events
    for j=1:org.numEvents
        
        % Assign the event number to a temporary variable
        event = ['Event_',num2str(j,'%04.0f')];
        
        % Find the first time code after the event in the full data file and
        % store the index
        x = find(raw.data.(name)(:,1) > org.eventTimeCodes(j),1);
        
        % Set the index range to capture
        start = x;
        stop = x + org.eventDuration(j)*org.sampleFreq;

        if stop > length(raw.data.(name))
            stop = length(raw.data.(name));
        end
        
        % Store the captured data into a temporary variable
        tempData = raw.data.(name)(start:stop,:);
        
        % Assign the temporary variable to a structured one organized by events
        rawSedentaryEvents.(event).(name) = tempData;
        clear x start stop tempData event
    end
    
    clear j name
end
clear i

%% Downsample large events
% this is done for program speed. Most large events will be removed as
% sleeping time, but regardless a sedentary bout recorded in hours does not
% need joint angles calculated at 20 Hz. This could be done for all
% sedentary bouts without significant loss if the computer running the
% program was unable to process at a reasonable speed.

% threshold data length
h = org.downsampleCutoff * 60 * 60 * org.sampleFreq;
org.downsampledEvents = [];

% Loop over the segments
for i=1:3
   
    % load the segment name as a temporary variable
    name = org.segmentNames{i};

    % Loop over the number of events
    for j=1:org.numEvents
        
        % Assign the event number to a temporary variable
        event = ['Event_',num2str(j,'%04.0f')];

        % Load temporary data
        tempData = rawSedentaryEvents.(event).(name);

        % Check to see if the event is longer than the desired cutoff
        if length(tempData) > h-1

            % If it is, downsample the data by the sample frequency,
            % lowering it to 1Hz
            tempDownsampled = downsample(tempData,org.sampleFreq);

            % Save the data
            rawSedentaryEvents.(event).(name) = tempDownsampled;

            % Create a list of which events have been downsampled
            org.downsampledEvents = [org.downsampledEvents event];
        end

        clear temp* event


    end
    clear j name

end
clear i h

%% Make the last event the same size

% Typically a final "event" is recorded when the APs are stationary prior
% to being downloaded. If this is the case the last event will be of
% different lengths in each AP as the thigh might have been turned off
% earlier than the others. Cropping the last event to be consistent with the
% shortest recording solves this. These events are ignored by the RA
% processing the output data.

% Find the last event
event = ['Event_',num2str(org.numEvents,'%04.0f')];

% Load the length of the last event from each AP
tempThigh = length(rawSedentaryEvents.(event).thigh);
tempTibia = length(rawSedentaryEvents.(event).tibia);
tempTorso = length(rawSedentaryEvents.(event).torso);

% Find the lowest value from the three
x = min([tempThigh tempTibia tempTorso]);

% Cut off values from the longer AP to match the shorter
newThigh = rawSedentaryEvents.(event).thigh(1:x,:);
newTibia = rawSedentaryEvents.(event).tibia(1:x,:);
newTorso = rawSedentaryEvents.(event).torso(1:x,:);

% Store the new data for the event
rawSedentaryEvents.(event).thigh = newThigh;
rawSedentaryEvents.(event).tibia = newTibia;
rawSedentaryEvents.(event).torso = newTorso;

clear temp* new* x event


%% Filter
%Low pass filter

Fs=org.sampleFreq; %sample freq
Fc=org.filterFreq; %filter freq
n=1; %filter order
Wn=Fc/(Fs/2);
[b,a] = butter(n,Wn, 'low'); 


% Loop over the number of events
for k=1:org.numEvents
    
    % Assign the event number to a temporary variable
    event = ['Event_',num2str(k,'%04.0f')];

    for i=1:3
       
        name = org.segmentNames{i};
    
        % Load the data from the respective event
        tempData = rawSedentaryEvents.(event).(name);
        
        % Initialize a temporary variable to store the filtered data
        tempFiltered = tempData;
        
        % Loop over the X Y and Z axis variables
        for j=1:3
        
            % Load the axis
            tempRowData = tempData(:,j+1);
        
            % Filter the data in that column
            tempRowFiltered = filtfilt(b,a,tempRowData);
        
            % Save the filtered data column
            tempFiltered(:,j+1) = tempRowFiltered;
        
            clear tempRow*
        end
        
        % Save the filtered data
        filteredSedentaryEvents.(event).(name) = tempFiltered;
        
        clear temp* j name

    end
    clear i event

end

clear k a b Fs Fc n Wn



%% Scale into Gs

% Scaling factors depending on activPAL version
if org.APversion == 3
    org.offset = 128;
    org.factor = 70;
elseif org.APversion == 4
    org.offset = 132;
    org.factor = 68;
else
    msg = 'Did not select a supported activPAL version';
    error(msg)
end

% Loop over the number of events
for k=1:org.numEvents
    
    % Assign the event number to a temporary variable
    event = ['Event_',num2str(k,'%04.0f')];
    
    for i=1:3
       
        name = org.segmentNames{i};
        
        % Load the data from the respective event
        tempData = filteredSedentaryEvents.(event).(name);
        
        tempScaled = tempData;
        
        for j=1:length(tempData)
            tempScaled(j,2) = (tempData(j,2)-org.offset)/org.factor;
            tempScaled(j,3) = (tempData(j,3)-org.offset)/org.factor;
            tempScaled(j,4) = (tempData(j,4)-org.offset)/org.factor;
        end
        
        clear j
        
        scaledSedentaryEvents.(event).(name) = tempScaled;
        
        clear temp* name
    end
    clear i event
end
clear k

%% Calculate pitch/Joint angle

f1 = waitbar(0, 'Starting Joint Angle calculation');
n = org.numEvents;

% Loop over the number of events
for k=1:org.numEvents
   
    
    % Assign the event number to a temporary variable
    event = ['Event_',num2str(k,'%04.0f')];

    % Dot product method to determine joint angles
    tempTibia = scaledSedentaryEvents.(event).tibia;
    tempThigh = scaledSedentaryEvents.(event).thigh;
    tempTorso = scaledSedentaryEvents.(event).torso;
    

    tempJoint = zeros(length(tempTibia),2);

    for j=1:length(tempTibia)

        

        a = tempTibia(j,2:4);
        b = tempThigh(j,2:4);
        c = tempTorso(j,2:4);

        tempJoint(j,1) = rad2deg(acos(dot(a/rssq(a),b/rssq(b))));
        tempJoint(j,2) = rad2deg(acos(dot(b/rssq(b),c/rssq(c))));
        
    end
    
    waitbar(k/n, f1, sprintf('Processing events: %d %%', floor(k/n*100)));

    dotProductJoints.(event) = tempJoint;
    clear j temp* event
%     close(f2)
end
clear k a b c

close(f1)
clear p n f1


%% Interpret joint values via various methods

% Pre-filter descriptive stats 

org.jointNames = {'knee','hip'};

% Loop over the number of events
for k=1:org.numEvents
    
    % Assign the event number to a temporary variable
    event = ['Event_',num2str(k,'%04.0f')];

    % loop over the calculated joints
    for j=1:2

        % Load a temporary variable with the joint name
        joint = org.jointNames{j};

        % Find the descriptive statistics for the joint angles
        jointAngle.(event).(joint).mean =...
            mean(dotProductJoints.(event)(:,j));
        jointAngle.(event).(joint).median =...
            median(dotProductJoints.(event)(:,j));
        jointAngle.(event).(joint).mode =...
            mode(dotProductJoints.(event)(:,j));
        jointAngle.(event).(joint).std =...
            std(dotProductJoints.(event)(:,j));


        % Find the proportion of each event within the angle bins
        proportionUnder15 = sum(dotProductJoints.(event)(:,j) < 15)...
            /length(dotProductJoints.(event)(:,j));
        proportionBtwn15_30 = sum(dotProductJoints.(event)(:,j) < 30 &...
            dotProductJoints.(event)(:,j) > 15)...
            /length(dotProductJoints.(event)(:,j));
        proportionBtwn30_45 = sum(dotProductJoints.(event)(:,j) < 45 &...
            dotProductJoints.(event)(:,j) > 30)...
            /length(dotProductJoints.(event)(:,j));
        proportionBtwn45_60 = sum(dotProductJoints.(event)(:,j) < 60 &...
            dotProductJoints.(event)(:,j) > 55)...
            /length(dotProductJoints.(event)(:,j));
        proportionBtwn60_75 = sum(dotProductJoints.(event)(:,j) < 75 &...
            dotProductJoints.(event)(:,j) > 60)...
            /length(dotProductJoints.(event)(:,j));
        proportionOver75 = sum(dotProductJoints.(event)(:,j) > 75)...
            /length(dotProductJoints.(event)(:,j));

        jointAngle.(event).(joint).secondsUnder15 =...
            org.eventDuration(k)*proportionUnder15;
        jointAngle.(event).(joint).secondsBtwn15_30 =...
            org.eventDuration(k)*proportionBtwn15_30;
        jointAngle.(event).(joint).secondsBtwn30_45 =...
            org.eventDuration(k)*proportionBtwn30_45;
        jointAngle.(event).(joint).secondsBtwn45_60 =...
            org.eventDuration(k)*proportionBtwn45_60;
        jointAngle.(event).(joint).secondsBtwn60_75 =...
            org.eventDuration(k)*proportionBtwn60_75;
        jointAngle.(event).(joint).secondsOver75 =...
            org.eventDuration(k)*proportionOver75;

        jointAngle.(event).(joint).proportionUnder15 = proportionUnder15;
        jointAngle.(event).(joint).proportionBtwn15_30 =...
            proportionBtwn15_30;
        jointAngle.(event).(joint).proportionBtwn30_45 =...
            proportionBtwn30_45;
        jointAngle.(event).(joint).proportionBtwn45_60 =...
            proportionBtwn45_60;
        jointAngle.(event).(joint).proportionBtwn60_75 =...
            proportionBtwn60_75;
        jointAngle.(event).(joint).proportionOver75 = proportionOver75;

        clear proportion* joint
    end

    clear event j
end
clear k

%% Output

% Creating a cell array that will store, in order: event number, event
% date&time, median knee angle, bin, and proportion of time spent in that
% bin, and the same three for hip angle
temp = cell(org.numEvents,17);

% Loop over the number of events
for i=1:org.numEvents
    
    % Assign the event number to a temporary variable
    event = ['Event_',num2str(i,'%04.0f')];

    % store the event number
    temp(i,1) = {i};
    
    % Find the initial time code for the event and convert to a date time
    % string
    temp(i,2) = {datestr(datetime(rawSedentaryEvents.(event).thigh(1,1),...
        "ConvertFrom",'excel'));};

    % Find the duration of the event and store in seconds
    s = seconds(org.eventDuration(i));
    % Reformat the duration into hours, minutes, seconds
    m = org.eventDuration(i)/60;
    s.Format = 'hh:mm:ss';
    % Store the duration
    temp(i,3) = {m};
    % temp(i,3) = {string(s, 'hh:mm:ss')};
    % Clear temp variable
    clear s m

    % Store the median knee joint angle for the event
    temp(i,4) = {jointAngle.(event).knee.median};

    % Time spent in < 15 bin
    s = seconds(jointAngle.(event).knee.secondsUnder15);
    m = jointAngle.(event).knee.secondsUnder15/60;
    s.Format = 'hh:mm:ss';
    temp(i,5) = {m};
    %temp(i,5) = {string(s, 'hh:mm:ss')};
    clear s m

    % Time spent in 15 - 30 bin
    s = seconds(jointAngle.(event).knee.secondsBtwn15_30);
    m = jointAngle.(event).knee.secondsBtwn15_30/60;
    s.Format = 'hh:mm:ss';
    temp(i,6) = {m};
    %temp(i,6) = {string(s, 'hh:mm:ss')};
    clear s m

    % Time spent in 30 - 45 bin
    s = seconds(jointAngle.(event).knee.secondsBtwn30_45);
    m = jointAngle.(event).knee.secondsBtwn30_45/60;
    s.Format = 'hh:mm:ss';
    temp(i,7) = {m};
    %temp(i,7) = {string(s, 'hh:mm:ss')};
    clear s m

    % Time spent in 45 - 60 bin
    s = seconds(jointAngle.(event).knee.secondsBtwn45_60);
    m = jointAngle.(event).knee.secondsBtwn45_60/60;
    s.Format = 'hh:mm:ss';
    temp(i,8) = {m};
    %temp(i,8) = {string(s, 'hh:mm:ss')};
    clear s m

    % Time spent in 60 - 75 bin
    s = seconds(jointAngle.(event).knee.secondsBtwn60_75);
    m = jointAngle.(event).knee.secondsBtwn60_75/60;
    s.Format = 'hh:mm:ss';
    temp(i,9) = {m};
    %temp(i,9) = {string(s, 'hh:mm:ss')};
    clear s m

    % Time spent in > 75 bin
    s = seconds(jointAngle.(event).knee.secondsOver75);
    m = jointAngle.(event).knee.secondsOver75/60;
    s.Format = 'hh:mm:ss';
    temp(i,10) = {m};
    %temp(i,10) = {string(s, 'hh:mm:ss')};
    clear s m

    temp(i,11) = {jointAngle.(event).hip.median};

    % Time spent in < 15 bin
    s = seconds(jointAngle.(event).hip.secondsUnder15);
    m = jointAngle.(event).hip.secondsUnder15/60;
    s.Format = 'hh:mm:ss';
    temp(i,12) = {m};
    %temp(i,12) = {string(s, 'hh:mm:ss')};
    clear s m

    % Time spent in 15 - 30 bin
    s = seconds(jointAngle.(event).hip.secondsBtwn15_30);
    m = jointAngle.(event).hip.secondsBtwn15_30/60;
    s.Format = 'hh:mm:ss';
    temp(i,13) = {m};
    %temp(i,13) = {string(s, 'hh:mm:ss')};
    clear s m

    % Time spent in 30 - 45 bin
    s = seconds(jointAngle.(event).hip.secondsBtwn30_45);
    m = jointAngle.(event).hip.secondsBtwn30_45/60;
    s.Format = 'hh:mm:ss';
    temp(i,14) = {m};
    %temp(i,14) = {string(s, 'hh:mm:ss')};
    clear s m

    % Time spent in 45 - 60 bin
    s = seconds(jointAngle.(event).hip.secondsBtwn45_60);
    m = jointAngle.(event).hip.secondsBtwn45_60/60;
    s.Format = 'hh:mm:ss';
    temp(i,15) = {m};
    %temp(i,15) = {string(s, 'hh:mm:ss')};
    clear s m

    % Time spent in 60 - 75 bin
    s = seconds(jointAngle.(event).hip.secondsBtwn60_75);
    m = jointAngle.(event).hip.secondsBtwn60_75/60;
    s.Format = 'hh:mm:ss';
    temp(i,16) = {m};
    %temp(i,16) = {string(s, 'hh:mm:ss')};
    clear s m*

    % Time spent in > 75 bin
    s = seconds(jointAngle.(event).hip.secondsOver75);
    m = jointAngle.(event).hip.secondsOver75/60;
    s.Format = 'hh:mm:ss';
    temp(i,17) = {m};
    %temp(i,17) = {string(s, 'hh:mm:ss')};
    clear s m*
    
    clear event

end
clear i

Output = table(temp(:,1),temp(:,2),temp(:,3),temp(:,4),temp(:,5),...
    temp(:,6),temp(:,7),temp(:,8),temp(:,9),temp(:,10),temp(:,11),...
    temp(:,12),temp(:,13),temp(:,14),temp(:,15),temp(:,16),temp(:,17),...
    'VariableNames',{'Event #','Date/Time','Duration (Minutes)',...
    'Median knee angle','Time knee spent < 15 (m)','Time knee spent in 15 - 30 (m)',...
    'Time knee spent in 30 - 45 (m)','Time knee spent in 45 - 60 (m)',...
    'Time knee spent in 60 - 75 (m)','Time knee spent over 75 (m)','Median hip angle',...
    'Time hip spent < 15 (m)','Time hip spent in 15 - 30 (m)',...
    'Time hip spent in 30 - 45 (m)','Time hip spent in 45 - 60 (m)',...
    'Time hip spent in 60 - 75 (m)','Time hip spent over 75 (m)'});

clear temp

%% output

% Reqest file save location and name from user, default based on origonal
% file
[file,path] = uiputfile({'*.xlsx'},'Save output as:','.xlsx');

% Write each output table as a sheet in the indicated file
writetable(Output,strcat(path,file),'Sheet','Sedentary Events');