%% SS.m
% Implements Seperate Systems Integration described by section 8.4.2 in the 
% report.
function [velocity, position, orientation] = ...
    SS(gyro, accel, mag, GPSPosition, ~, ...
    initialPosition, initialVelocity, initialOrientation, delta,...
    rateRatio, globalG, globalM)
    
    % GPS samples between each reset
    resetPeriod = 100;
    
    % Initialise   
    gpsFS = 1/delta/rateRatio;
    currVelocity = initialVelocity;
    currOrientation = initialOrientation;
    currPosition = initialPosition;
    velocity = [];
    position = [];
    orientation = [];    

    for i = 1 : length(GPSPosition) / resetPeriod - 1
        s = i - 1;
        blockRows = (rateRatio*resetPeriod*s + 1):...
            (rateRatio*resetPeriod*s+rateRatio*resetPeriod);
        gyroBlock = gyro(blockRows,:);
        accelBlock = accel(blockRows,:);
        magBlock = mag(blockRows,:);
        gpsBlock = GPSPosition(s*resetPeriod + 1:s*resetPeriod + resetPeriod,:);
        [velocityTemp, positionTemp, orientationTemp] = ...
            simple(gyroBlock, accelBlock, 0, 0, 0, ...
            currPosition, currVelocity, currOrientation, delta, 0);  

        % Linear fit to gps positions to find velocity vector
        fitX = polyfit(linspace(0,resetPeriod/gpsFS,resetPeriod),gpsBlock(:,1),1);  
        fitY = polyfit(linspace(0,resetPeriod/gpsFS,resetPeriod),gpsBlock(:,2),1);  
        fitZ = polyfit(linspace(0,resetPeriod/gpsFS,resetPeriod),gpsBlock(:,3),1);

        % Reset current position, velocity and orientation.
        currPosition = gpsBlock(end,:);
        currVelocity = [fitX(1),fitY(1),fitZ(1)];

        % Orientation Estimation
        q = ...
        quaternionAttitudeEstimator(globalG,mean(magBlock)',globalG,globalM,1,1);

        % Orientation Reset
        currOrientation = q;    

        % Output Position Update
        velocity = [velocity;velocityTemp];
        position = [position;positionTemp];
        orientation = [orientation;orientationTemp];
    end    
end

