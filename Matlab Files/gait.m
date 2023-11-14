% Importing (reading) data
AccData = readtable('pip_gait.csv');
%AccData = AccData(AccData(:,4) == "s11p1c1_B1.trc", :);
% Assuming AccData is a matrix with numeric values in the fourth column
% Assuming AccData is a table with a string column named 'Column4'
targetString = 's11p1c1_B1.trc';

% Filter the rows based on the fourth column value
AccData = AccData(strcmp(AccData.subject, targetString), :);


%AccData = AccData(strcmp(AccData(:,4), 's11p1c1_B1.trc'), :);

%plot(AccData(:,3));
%AccData = AccData(197:597,:);
% Integration and peak detection 
fs = 60;
%av = AccData(:,2); 
av = AccData{:,2};
%d = class(av);

Integratedav = cumtrapz(av); % Numerical integration of av
s1 = 1*(cwt(Integratedav, 14, 'gaus1', 1/fs)); % using Gaussian CWT at scale 10 AND 14 FOR C34 AND C57 AND P445 AND P466 AND P487 AND P503
[Peaks, Locations1] = findpeaks(s1); % Locates the minimum Peaks (values) and locations (time/samples)
IC = Locations1;
s1 = -1*s1,
s2 = cwt(s1, 14, 'gaus1', 1/fs);
[Peaks2, Locations2] = findpeaks(s2);
FC = Locations2;

axvec = 1 : length(av);
minav = min(av); mins1 = min(s1); mins2 = min(s2);
maxav = max(av); maxs1 = max(s1); maxs2 = max(s2); 
xs = (av(:) - minav) ./ (maxav - minav);  
ys = (s1(:) - mins1) ./ (maxs1 - mins1);
zs = (s2(:) - mins2) ./ (maxs2 - mins2);
%plot(axvec, xs, axvec, ys, axvec, zs);
%legend('Av','S1','S2'); 
%xlabel('Samples'); 

for i = 2:length(FC) % The for loop must be iterate a predefined
 StepTime(i-1) = IC(i) - IC(i-1);
 StanceTime(i-1) = FC(i) - IC(i); 
end

for i = 2:length(FC)
StrideTime(i-1) = IC(i + 1) - IC(i-1);
end
StanceTime1 = StanceTime;
%StanceTime1(end) = [];
StrideTime1 = StrideTime;
%StrideTime1(end) = [];
StrideTime1 = StrideTime1(1:end-1);

SwingTime = StrideTime1 - StanceTime1;

StepTime = StepTime.'
StanceTime = StanceTime.'
StrideTime = StrideTime.'
SwingTime = SwingTime.'

outputData = [AccData.subject, StepTime, StanceTime, StrideTime, SwingTime];
outputFilename = 's11p1c1_B1.trc.csv';
csvwrite(outputFilename, outputData);
disp('Saved.');

% Variability and asymmetry calculations
%StepTimeLeft = StepTime(1:2:end,:); % from 1st data point to end
%StepTimeLeft1 = StepTimeLeft.'
%StepTimeRight = StepTime(2:2:end,:); % from 2nd data point to end
%StepTimeRight1 = StepTimeRight.'
% Step time 
%StepTimeVrl = sqrt((var(StepTimeLeft1)+var(StepTimeRight1))/2); % Equation 7a
%StepTimeV = std(StepTime); % Equation 7b
%StepTimeA = abs(mean(StepTimeLeft1) - mean(StepTimeRight1)); % Equation 8
% Stance time 
%StanceTimeLeft = StanceTime(1:2:end,:); % from 1st data point to end
%StanceTimeLeft1 = StanceTimeLeft.'
%StanceTimeRight = StanceTime(2:2:end,:); % from 2nd data point to end
%StanceTimeRight1 = StanceTimeRight.'
%StanceTimeVrl = sqrt((var(StanceTimeLeft1)+var(StanceTimeRight1))/2); % Equation 7a
%StanceTimeV = std(StanceTime); % Equation 7b
%StanceTimeA = abs(mean(StanceTimeLeft1) - mean(StanceTimeRight1)); % Equation 8
% Stride time 
%StrideTimeLeft = StrideTime(1:2:end,:); % from 1st data point to end
%StrideTimeLeft1 = StrideTimeLeft.'
%StrideTimeRight = StrideTime(2:2:end,:); % from 2nd data point to end
%StrideTimeRight1 = StrideTimeRight.'
%StrideTimeVrl = sqrt((var(StrideTimeLeft1)+var(StrideTimeRight1))/2); % Equation 7a
%StrideTimeV = std(StrideTime); % Equation 7b
%StrideTimeA = abs(mean(StrideTimeLeft1) - mean(StrideTimeRight1)); % Equation 8
% Swing time 
%SwingTimeLeft = SwingTime(1:2:end,:); % from 1st data point to end
%SwingTimeLeft1 = SwingTimeLeft.'
%SwingTimeRight = SwingTime(2:2:end,:); % from 2nd data point to end
%SwingTimeRight1 = SwingTimeRight.'
%SwingTimeVrl = sqrt((var(SwingTimeLeft1)+var(SwingTimeRight1))/2); % Equation 7a
%SwingTimeV = std(SwingTime); % Equation 7b
%SwingTimeA = abs(mean(SwingTimeLeft1) - mean(SwingTimeRight1)); % Equation 8

%aaa = [StepTime StanceTime];
%bbb = [StrideTime SwingTime];
%ccc = [StepTimeVrl StepTimeV	StepTimeA	StanceTimeVrl	StanceTimeV	StanceTimeA	StrideTimeVrl	StrideTimeV	StrideTimeA	SwingTimeVrl	SwingTimeV	SwingTimeA];