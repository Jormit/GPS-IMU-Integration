%% allan.m
% This code implements the Allan Deviation method of characterising IMU
% noise characteristics as described by section 7.2.2 in the report. By
% default this is using data from the file "ADIS16507 Gyro 8 minutes.csv".
% This code is based upon this example:
% https://au.mathworks.com/help/fusion/ug/inertial-sensor-noise-analysis-using-allan-variance.html

%% Setup
fs = 1000;
data = readmatrix("ADIS16507 Gyro 8 minutes.csv");
omega = data - mean(data);

%% Allan Deviation Plot
[avar, tau] = sigTau(omega, fs);
adev = sqrt(avar);
figure
loglog(tau, adev, 'LineWidth',2); 
title('Allan Deviation')
xlabel('\tau');
ylabel('\sigma(\tau)')
grid on
hold on

%% Angle Random Walk Determination
slope = -0.5;
[index, logadev, logtau] = findAllanSlope(adev, tau, slope);
b = logadev(index) - slope*logtau(index);
logAWR = slope*log10(1) + b;
awr = 10^logAWR;
lineAWR = awr ./ sqrt(tau);
loglog(tau,lineAWR,"--")
loglog(1,awr,"o")

%% Rate Random Walk Determination
slope = 0.5;
[index, logadev, logtau] = findAllanSlope(adev, tau, slope);
b = logadev(index) - slope*logtau(index);
logRRW = slope*log10(3) + b;
rrw = 10^logRRW;
lineRRW = rrw .* sqrt(tau/3);
loglog(tau,lineRRW,"--")
loglog(3,rrw,"o")

%% Bias Instability Determination
slope = 0;
[index, logadev, logtau] = findAllanSlope(adev, tau, slope);
b = logadev(index) - slope*logtau(index);
scfBI = sqrt(2*log(2)/pi);
logBI = b - log10(scfBI);
bi = 10^logBI;
lineBI = bi * scfBI * ones(size(tau));
loglog(tau,lineBI,"--")

legend('Allan Deviation',...
       'Angle Random Walk Slope',...
       'Angle Random Walk \tau = 1 Intercept',...
       'Rate Random Walk Slope',...
       'Rate Random Walk \tau = 3 Intercept',...
       'Bias Instability Slope')

hold off;

%% Model Comparison

% Simulate using identified parameters
IMU = imuSensor('accel-gyro'); 
IMU.Gyroscope = gyroparams( ...
        'NoiseDensity',awr, ...
        'BiasInstability',bi, ...
        'RandomWalk',rrw);  
IMU.SampleRate = fs;
[~,gyroReading] = IMU(zeros(length(data),3),zeros(length(data),3));

[avar1, tau1] = sigTau(omega, fs);
[avar2, tau2] = sigTau(gyroReading(:,1), fs);

% Compare Allan Plots
figure
loglog(tau1, sqrt(avar1), tau2, sqrt(avar2));
title('Allan Deviation Comparison')
xlabel('\tau');
ylabel('\sigma(\tau)')
grid on
hold on
legend('Measured', 'Simulated');

% Compare Time Domain Plots
figure
subplot(1,2,1)
plot(omega);
title('Measured')
xlabel('Samples');
ylabel('Angular Rate (\circ/sec)')
subplot(1,2,2)
plot(gyroReading(:,1));
title('Simulated')
xlabel('Samples');
ylabel('Angular Rate (\circ/sec)')

% Compare Against Datasheet
biasInstability = [bi*3600; 7.5];
angularRandomWalk = [awr*60; 0.32];
labels = ["Calculated"; "Datasheet"];
table(labels, biasInstability, angularRandomWalk)

%% Allan Variance Sigma Tau
function [sigma, tau] = sigTau(omega, fs)
    maxNumM = 100;
    L = size(omega, 1);
    maxM = 2.^floor(log2(L/2));
    m = logspace(log10(1), log10(maxM), maxNumM).';
    m = ceil(m);
    m = unique(m);
    [sigma, tau] = allanvar(omega, m, fs);
end

%% Returns the index of allandev slope closest to specified slope.
function [index, logadev, logtau] = findAllanSlope(adev, tau, slope)
    logtau = log10(tau);
    logadev = log10(adev);
    dlogadev = diff(logadev) ./ diff(logtau);
    [~, index] = min(abs(dlogadev - slope));
end   