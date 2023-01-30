%% EKF.m
% Implements Extended Kalman Filter with bias estimation described by 
% section 8.4.5 and 7.3.5.9 in the report.
function [velocity, position, orientation] = ...
    EKFBias(gyro, accel, mag, GPSPosition, GPSVelocity, ...
    initialPosition, initialVelocity, initialOrientation, delta,...
    rateRatio, globalG, globalM)

    % Field Parameters
    g = 9.81;

    % Sensor Parameters
    GPSPosVar = 150/rateRatio;
    GPSVelVar = 200/rateRatio;
    AccVar = 2.4314e-04;
    RateVar = 4.9013e-06;
    MagVar = 3.4138e-05;
    BiasVar = 0;

    % Initial State Vector
    X = [initialPosition';initialVelocity';initialOrientation';zeros(6,1)];

    % Initial State Covariance
    P = eye(16);

    for i = 1 : length(accel)
        % IMU Readings
        w = gyro(i,:)';
        a = accel(i,:)';
        m = mag(i,:)';

        % State Prediction
        omega = Omega(w-X(14:16));
        q = X(7:10);
        Rq = quat2rotm(X(7:10)');
        dvdq = dxdq((a - X(11:13)),X(7:10));
        F = eye(16) + [zeros(3)  ,eye(3)    ,zeros(3,4), zeros(3)  , zeros(3);
                       zeros(3)  ,zeros(3)  ,dvdq      , -Rq       , zeros(3);
                       zeros(4,3),zeros(4,3),0.5*omega , zeros(4,3), -xi(q)/2;
                       zeros(3)  ,zeros(3)  ,zeros(3,4), zeros(3)  , zeros(3);
                       zeros(3)  ,zeros(3)  ,zeros(3,4), zeros(3)  , zeros(3)] * delta;    
        Xp = [X(1:3) + X(4:6) * delta;
              X(4:6) + (Rq * (a - X(11:13)) - [0,0,g]') * delta;
              X(7:10) + 0.5 * omega * delta * X(7:10);
              X(11:13);
              X(14:16)];

        % Covariance prediction
        W = 0.5 * delta * xi(q);
        posCov  = AccVar * delta^4 * eye(3);
        velCov  = AccVar * delta^2 * eye(3);
        RateCov = RateVar * (W) * (W');
        BiasCov = BiasVar * eye(3);
        Q = [posCov    ,zeros(3)  ,zeros(3,4),zeros(3)  ,zeros(3)  ;
             zeros(3)  ,velCov    ,zeros(3,4),zeros(3)  ,zeros(3)  ;
             zeros(4,3),zeros(4,3),RateCov   ,zeros(4,3),zeros(4,3);
             zeros(3)  ,zeros(3)  ,zeros(3,4),BiasCov   ,zeros(3)  ;
             zeros(3)  ,zeros(3)  ,zeros(3,4),zeros(3)  ,BiasCov  ];     
        Pp = F * P * F' + Q;

        % Measurement Matrix and Covariance
        qp = Xp(7:10);
        qpc = quatconj(qp')';
        gMeas = normalise(a);
        gMeas = globalG;
        mMeas = normalise(m);
        dgdq = dxdq(globalG,qpc);    
        dmdq = dxdq(globalM,qpc);
        C = quat2rotm(qp')';

        GPSReady = ~mod(i-1, rateRatio);
        z = [GPSPosition(floor((i-1)/rateRatio) + 1, :)';...
             GPSVelocity(floor((i-1)/rateRatio) + 1, :)';...
             gMeas;mMeas];
        h = [Xp(1:6); C*globalG; C*globalM];
        H = [eye(3)  ,zeros(3),zeros(3,4),zeros(3),zeros(3);
             zeros(3),eye(3)  ,zeros(3,4),zeros(3),zeros(3);
             zeros(3),zeros(3),dgdq      ,zeros(3),zeros(3);
             zeros(3),zeros(3),dmdq      ,zeros(3),zeros(3)];
        R = [GPSPosVar*eye(3),zeros(3)        ,zeros(3)     ,zeros(3)     ;
             zeros(3)        ,GPSVelVar*eye(3),zeros(3)     ,zeros(3)     ;
             zeros(3)        ,zeros(3)        ,AccVar*eye(3),zeros(3)     ;
             zeros(3)        ,zeros(3)        ,zeros(3)     ,MagVar*eye(3)];       

        % Kalman Gain
        K = Pp * H' / (H * Pp * H' + R);
        
        % Ignore GPS readings if not ready.
        if (~GPSReady)
            K(:,1:3)=0;
            K(:,4:6)=0;
        end 

        % State Update
        X = Xp + K * (z - h);

        % State Covariance Update
        P = (eye(16) - K * H) * Pp;

        % Normalise Quaternion
        X(7:10) = normalise(X(7:10));

        % Save state vector.
        position(i, :) = X(1:3)';
        velocity(i, :) = X(4:6)';
        orientation(i, :) = X(7:10)';
        accelBias(i, :) = X(11:13)';
        gyroBias(i, :) = X(14:16)';
        accelbiasConfidence(i, :) = sqrt([P(11,11),P(12,12),P(13,13)]);
        gyrobiasConfidence(i, :) = sqrt([P(14,14),P(15,15),P(16,16)]);
    end
    
    % Accel Bias confidence intervals
    figure
    subplot(3,2,1)
    x = linspace(0,200,length(accel));                     
    y = accelBias(:,1)';
    c = 3*accelbiasConfidence(:,1)';
    xconf = [x x(end:-1:1)] ;         
    yconf = [y+c y(end:-1:1)-c(end:-1:1)];
    p = fill(xconf,yconf,'red');
    p.FaceColor = [1 0.8 0.8];      
    p.EdgeColor = 'none';
    hold on
    plot(x,y,'r')
    xlim([1,40]);
    hold off
    legend('3\sigma Bound','Estimated Bias')
    xlabel('Time (s)')
    title('X Accelerometer Bias Estimate')
    ylabel('Bias (m/s^2)')

    subplot(3,2,3)
    x = linspace(0,200,length(accel));                     
    y = accelBias(:,2)';
    c = 3*accelbiasConfidence(:,2)';
    xconf = [x x(end:-1:1)] ;         
    yconf = [y+c y(end:-1:1)-c(end:-1:1)];
    p = fill(xconf,yconf,'red');
    p.FaceColor = [1 0.8 0.8];      
    p.EdgeColor = 'none';
    hold on
    plot(x,y,'r')
    xlim([1,40]);
    hold off
    legend('3\sigma Bound','Estimated Bias')
    xlabel('Time (s)')
    title('Y Accelerometer Bias Estimate')
    ylabel('Bias (m/s^2)')

    subplot(3,2,5)
    x = linspace(0,200,length(accel));                     
    y = accelBias(:,3)';
    c = 3*accelbiasConfidence(:,3)';
    xconf = [x x(end:-1:1)] ;         
    yconf = [y+c y(end:-1:1)-c(end:-1:1)];
    p = fill(xconf,yconf,'red');
    p.FaceColor = [1 0.8 0.8];      
    p.EdgeColor = 'none';
    hold on
    plot(x,y,'r')
    xlim([1,40]);
    hold off
    legend('3\sigma Bound','Estimated Bias')
    xlabel('Time (s)')
    title('Z Accelerometer Bias Estimate')
    ylabel('Bias (m/s^2)')

    % Gyro Bias confidence intervals
    subplot(3,2,2)
    x = linspace(0,200,length(accel));                     
    y = gyroBias(:,1)';
    c = 3*gyrobiasConfidence(:,1)';
    xconf = [x x(end:-1:1)] ;         
    yconf = [y+c y(end:-1:1)-c(end:-1:1)];
    p = fill(xconf,yconf,'red');
    p.FaceColor = [1 0.8 0.8];      
    p.EdgeColor = 'none';
    hold on
    plot(x,y,'r')
    xlim([1,40]);
    hold off
    legend('3\sigma Bound','Estimated Bias')
    xlabel('Time (s)')
    title('X Gyroscope Bias Estimate')
    ylabel('Bias (rads/sec)')

    subplot(3,2,4)
    x = linspace(0,200,length(accel));                     
    y = gyroBias(:,2)';
    c = 3*gyrobiasConfidence(:,2)';
    xconf = [x x(end:-1:1)] ;         
    yconf = [y+c y(end:-1:1)-c(end:-1:1)];
    p = fill(xconf,yconf,'red');
    p.FaceColor = [1 0.8 0.8];      
    p.EdgeColor = 'none';
    hold on
    plot(x,y,'r')
    xlim([1,40]);
    hold off
    legend('3\sigma Bound','Estimated Bias')
    xlabel('Time (s)')
    title('Y Gyroscope Bias Estimate')
    ylabel('Bias (rads/sec)')
    
    subplot(3,2,6)
    x = linspace(0,200,length(accel));                     
    y = gyroBias(:,3)';
    c = 3*gyrobiasConfidence(:,3)';
    xconf = [x x(end:-1:1)] ;         
    yconf = [y+c y(end:-1:1)-c(end:-1:1)];
    p = fill(xconf,yconf,'red');
    p.FaceColor = [1 0.8 0.8];      
    p.EdgeColor = 'none';
    hold on
    plot(x,y,'r')
    xlim([1,40]);
    hold off
    legend('3\sigma Bound','Estimated Bias')
    xlabel('Time (s)')
    title('Z Gyroscope Bias Estimate')
    ylabel('Bias (rads/sec)')
end