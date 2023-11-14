% Importing (reading) data
AccData = xlsread('TAS1H06190087(2021-02-15)RAW.csv');
%plot(AccData(:,1));
% Plotting and manual segmentation
plot(AccData); % plots tri-axial accelerometer data
legend('Vertical','Medio','Antero'); % attaches a legend to the plot figure/window
xlabel('Samples'); ylabel('Gravity (g)'); % labels x and y
StartStop = ginput(2); % saves 2 x and y coordinates from user input clicking on plot
StartStop(:,1);
%AccData = AccData(300302:857800,:);

%AccData = AccData(444000:474000,:);

%AccData = AccData(684000:714000,:);

AccData = AccData(690000:714000,:);


% Filtering (It should be enough to do the process for only 1 axis(?)
fc = 20; % with a cut off frequency of 15–20Hz
fs = 100; % Define sample frequency as 100Hz
Wn = fc/(fs/2);
[B,A] = butter(4, Wn);

AccDataV = AccData(:,2); 
AccDataFilt1 = filtfilt(B, A, AccDataV);
%freqz(B,A);

AccDataM = AccData(:,1); 
AccDataFilt2 = filtfilt(B, A, AccDataM);

AccDataA = AccData(:,3); 
AccDataFilt3 = filtfilt(B, A, AccDataA);

AccData = [AccDataFilt1 AccDataFilt2 AccDataFilt3];
%plot(AccData);


% Acceleration correction to horizontal-vertical frame
av = AccData(:,1);
aa = AccData(:,3); 
am = AccData(:,2); 
aaMean = mean(aa);
avMean = mean(av); % Mean of av as best estimate of sin(?v)
amMean = mean(am); % Mean of am best estimate of sin(?m)
aA = aa*cos(asin(aaMean)) - av.*aaMean; % Equation 1
avv = aa*aaMean + asin(avMean).*cos(asin(aaMean)); % Equation 2 
aM = am*cos(asin(amMean)) - avv.*amMean; % Equation 3
aV = am.*amMean + avv.*cos(asin(amMean)) - 1; % Equation 4
AccData = [aV aM aA];
%plot(AccData);


% Integration and peak detection 
av = AccData(:,1); 
Integratedav = cumtrapz(av); % Numerical integration of av
s1 = -1*(cwt(Integratedav, 46, 'gaus1', 1/fs)); % using Gaussian CWT at scale 10
[Peaks, Locations1] = findpeaks(s1); % Locates the minimum Peaks (values) and locations (time/samples)
IC = Locations1;
s2 = cwt(s1, 55000, 'gaus1', 1/fs);
[Peaks2, Locations2] = findpeaks(s2);
FC = Locations2;
test = [IC FC];
plot(test);
%tt = s1.';
%tt = tt(1:557400,:);
%plot(tt);

% Calculating temporal gait characteristics from the IC/FC events
%for i = 3:length(IC) % The for loop must be iterate a predefined
 %StepTime(i-2) = IC(i-1) - IC(i-2);
 %StanceTime(i-2) = FC(i-1) - IC(i-2);
 %StrideTime(i-2) = IC(i) - IC(i-2);
%end
%SwingTime = StrideTime - StanceTime;

for i = 1:length(IC)-2 % The for loop must be iterate a predefined
 StepTime(i) = IC(i + 1) - IC(i);
 StanceTime(i) = FC(i + 1) - IC(i); %It was FC(i + 1) - IC(i);
 StrideTime(i) = IC(i + 2) - IC(i);
end
StepTime = StepTime(:,1:1434);
StanceTime = StanceTime(:,1:1434);
StrideTime = StrideTime(:,1:1434);
SwingTime = StrideTime - StanceTime;
plot(StepTime);
plot(StanceTime);
plot(StrideTime);
plot(SwingTime);

% Integration and spatio-temporal estimations
WearableHeight = 96;
hvel = cumtrapz(av); % Integrate to estimate velocity,
h = cumtrapz(hvel); % further integration will derive position
%StepLength = 2*(sqrt(abs(2*(WearableHeight)*h - h.^2))); % Inverted pendulum, Equation 5
StepLength = 22*(WearableHeight)*h - h.^2;
StepVelocity = StepLength/StepTime;
len = abs(2*(WearableHeight)*h- h.^2);

% Variability and asymmetry calculations
StepTimeLeft = StepTime(1:2:end,:); % from 1st data point to end
StepTimeRight = StepTime(2:2:end,:); % from 2nd data point to end
% Step time variability & asymmetry:
StepTimeV = sqrt((var(StepTimeLeft)+var(StepTimeRight))/2); % Equation 7a
% or
StepTimeV = std(StepTime); % Equation 7b
% and
StepTimeA = abs(mean(StepTimeLeft) - mean(StepTimeRight)); % Equation 8
