data = xlsread('P432.csv');
%data = xlsread('NAF_t05_untrimed.csv');
%data = data(121:402,:);

dt = 1/60; 
dhx = diff(data(:,1))
dhy = diff(data(:,2))
dhz = diff(data(:,3))

vx = dhx./dt
vy = dhy./dt
vz = dhz./dt
v = [vx vy vz];

dvx = diff(v(:,1))
dvy = diff(v(:,2))
dvz = diff(v(:,3))

ax = dvx./dt
ay = dvy./dt
az = dvz./dt

% filtering acceleration 
fc = 15; % with a cut off frequency of 15–20Hz
fs = 60; % Define sample frequency as 100Hz
Wn = fc/(fs/2);
[B,A] = butter(4, Wn);
ax1 = filtfilt(B, A, ax);
ay2 = filtfilt(B, A, ay); 
az3 = filtfilt(B, A, az);
a = [ax1 ay2 az3];
%plot(az3);
a = a(500:4100,:);

%csvwrite('NAF_t05_filtered.csv',a)
csvwrite('P432_filtered.csv',a)