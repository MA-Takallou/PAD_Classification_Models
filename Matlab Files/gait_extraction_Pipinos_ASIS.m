AccData = readtable('pip_gait.csv');

% Assuming the subject identifier is in the first column of the data
subjectIDs = unique(cellstr(AccData(:,4))); % Get unique subject IDs

% Create empty arrays to store the results for each subject
subjectStepTime = cell(size(subjectIDs));
subjectStanceTime = cell(size(subjectIDs));
subjectStrideTime = cell(size(subjectIDs));
subjectSwingTime = cell(size(subjectIDs));

for subjectIndex = 1:length(subjectIDs)
    subjectData = AccData(strcmpi(AccData(:,4), subjectIDs{subjectIndex}), :);

 
    fs = 60;
    Integratedav = cumtrapz(subjectData); % Numerical integration of av
    s1 = 1*(cwt(Integratedav, 10, 'gaus1', 1/fs)); % using Gaussian CWT at scale 10 AND 14 FOR C34 AND C57 AND P445 AND P466 AND P487 AND P503
    [Peaks, Locations1] = findpeaks(s1); % Locates the minimum Peaks (values) and locations (time/samples)
    IC = Locations1;
    s1 = -1*s1,
    s2 = cwt(s1, 10, 'gaus1', 1/fs);
    [Peaks2, Locations2] = findpeaks(s2);
    FC = Locations2;

    axvec = 1 : length(subjectData);
    minav = min(av); mins1 = min(s1); mins2 = min(s2);
    maxav = max(av); maxs1 = max(s1); maxs2 = max(s2); 
    xs = (av(:) - minav) ./ (maxav - minav);  
    ys = (s1(:) - mins1) ./ (maxs1 - mins1);
    zs = (s2(:) - mins2) ./ (maxs2 - mins2); 

    for i = 2:length(FC) % The for loop must be iterate a predefined
     StepTime(i-1) = IC(i) - IC(i-1);
     StanceTime(i-1) = FC(i) - IC(i-1); 
    end

    for i = 2:length(FC)-1
    StrideTime(i-1) = IC(i + 1) - IC(i-1);
    end
    StanceTime1 = StanceTime;
%StanceTime1(end) = [];
    StrideTime1 = StrideTime;
%StrideTime1(end) = [];
    SwingTime = StrideTime1 - StanceTime1;
    % Calculate gait parameters
    % ...
    % Your existing code to calculate StepTime, StanceTime, StrideTime, SwingTime goes here
    % ...

    % Append the results for the current subject to the corresponding arrays
    subjectStepTime = [subjectStepTime; StepTime.'];
    subjectStanceTime = [subjectStanceTime; StanceTime.'];
    subjectStrideTime = [subjectStrideTime; StrideTime.'];
    subjectSwingTime = [subjectSwingTime; SwingTime.'];
end

% Save the results for all subjects in a single CSV file
outputData = [subjectIDs, subjectStepTime, subjectStanceTime, subjectStrideTime, subjectSwingTime];
outputFilename = 'gait_results.csv';
csvwrite(outputFilename, outputData);
disp('Results saved successfully.');
