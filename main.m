%% main.m
%  This file is the driver code to test each of the various IMU/GPS
%  integration algorithms. It uses the trajectory defined in
%  generateTrajectory.m and generates simulated GPS and IMU measurements
%  with parameters defined in gpsModel.m and imuModel.m. By default these
%  are set to be representative of the ADIS16705 IMU and a GPS receiver
%  accurate to 10m. 

%% Create Trajectory
imuFS = 100;
delta = 1/imuFS;
[position,orientation,velocity,acceleration,angularVelocity, referenceLocation] = ...
    generateTrajectory(imuFS);

%% Simulate GPS Sensor
gpsFS = 10;
windowSize = 20;
GPSPositionLLA = ... 
    gpsModel(gpsFS, imuFS, position, velocity, referenceLocation);
GPSPosition = geo2enu(GPSPositionLLA,referenceLocation,wgs84Ellipsoid);
GPSVelocity = [0,0,0; diff(movmean(GPSPosition, windowSize))*gpsFS];

%% Simulate Inertial Sensor
globalG = [0,0,1]';
globalM = [0,1,0]';
[accelReading,gyroReading,magReading] = imuModel(acceleration,angularVelocity,orientation,imuFS);

%% Pass to integration algorithm
rateRatio = imuFS/gpsFS;

% Change this function to either EKFBias, EKF, UKF, SS and simple
[fusedVelocity, fusedPosition, fusedOrientation] = ... 
    EKF(gyroReading, accelReading, magReading, GPSPosition, GPSVelocity, ...
    position(1,:), velocity(1,:), compact(orientation(1)), delta,...
    rateRatio, globalG, globalM);

%% Plot Results
%% GPS vs GPS/IMU 2D Position
figure
plot(GPSPosition(:,1), GPSPosition(:,2), 'r.');
hold on
plot(fusedPosition(:,1), fusedPosition(:,2), 'b.');
plot(position(:,1), position(:,2), 'g')
legend('GPS Position', 'GPS/IMU Position', 'True Position')
ylabel('Y Position (m)');
xlabel('X Position (m)');
daspect([1,1,1])

%% True vs GPS/IMU Orientation
figure
subplot(2,1,1)
plot(rad2deg(quat2eul(compact(orientation))))
legend('z','y','x')
xlabel('Time (s)');
ylabel('Angle (\circ)');
title('True Orientation')
subplot(2,1,2)
plot(rad2deg(quat2eul(fusedOrientation)))
legend('z','y','x')
xlabel('Time (s)');
ylabel('Angle (\circ)');
title('GPS/IMU Orientation')

%% GPS vs GPS/IMU Position Error
figure
subplot(3,1,1);
hold on
title('GPS/IMU vs GPS Position Error (X-axis)')
plot(GPSPosition(:,1) - position(1:rateRatio:end,1))
plot(fusedPosition(1:rateRatio:end,1) - position(1:rateRatio:length(fusedPosition),1))
legend('GPS error', 'EKF error')
subplot(3,1,2);
hold on
title('GPS/IMU vs GPS Position Error (Y-axis)')
plot(GPSPosition(:,2) - position(1:rateRatio:end,2))
plot(fusedPosition(1:rateRatio:end,2) - position(1:rateRatio:length(fusedPosition),2))
legend('GPS error', 'EKF error')
subplot(3,1,3);
hold on
title('GPS/IMU vs GPS Position Error (Z-axis)')
plot(GPSPosition(:,3) - position(1:rateRatio:end,3))
plot(fusedPosition(1:rateRatio:end,3) - position(1:rateRatio:length(fusedPosition),3))
legend('GPS error', 'EKF error')

%% GPS vs GPS/IMU Velocity Error
figure
subplot(3,1,1);
hold on
title('EKF vs GPS Velocity Error (X-axis)')
plot(GPSVelocity(:,1) - velocity(1:rateRatio:end,1))
plot(fusedVelocity(1:rateRatio:end,1) - velocity(1:rateRatio:length(fusedVelocity),1))
legend('GPS error', 'EKF error')
subplot(3,1,2);
hold on
title('EKF vs GPS Velocity Error (Y-axis)')
plot(GPSVelocity(:,2) - velocity(1:rateRatio:end,2))
plot(fusedVelocity(1:rateRatio:end,2) - velocity(1:rateRatio:length(fusedVelocity),2))
legend('GPS error', 'EKF error')
subplot(3,1,3);
hold on
title('EKF vs GPS Velocity Error (Z-axis)')
plot(GPSVelocity(:,3) - velocity(1:rateRatio:end,3))
plot(fusedVelocity(1:rateRatio:end,3) - velocity(1:rateRatio:length(fusedVelocity),3))
legend('GPS error', 'EKF error')

%% Raw IMU measurements
figure
subplot(2,3,1)
plot(accelReading(:,1))
title('X-Axis Accleration')
xlabel('Sample Number')
ylabel('Acceleration (m/s^2)')
subplot(2,3,2)
plot(accelReading(:,2))
title('Y-Axis Accleration')
xlabel('Sample Number')
ylabel('Acceleration (m/s^2)')
subplot(2,3,3)
plot(accelReading(:,3))
title('Z-Axis Accleration')
xlabel('Sample Number')
ylabel('Acceleration (m/s^2)')
subplot(2,3,4)
plot(gyroReading(:,1))
title('X-Axis Angular Rate')
xlabel('Sample Number')
ylabel('Angular Rate (rads/sec)')
subplot(2,3,5)
plot(gyroReading(:,2))
title('Y-Axis Angular Rate')
xlabel('Sample Number')
ylabel('Angular Rate (rads/sec)')
subplot(2,3,6)
plot(gyroReading(:,3))
title('Z-Axis Angular Rate')
xlabel('Sample Number')
ylabel('Angular Rate (rads/sec)')
