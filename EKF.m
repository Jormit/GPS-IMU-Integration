%% EKF.m
% Implements Extended Kalman Filter described by section 8.4.3 and 7.3.5.9
% in the report.
function [velocity, position, orientation] = ...
    EKF(gyro, accel, mag, GPSPosition, GPSVelocity, ...
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
    gravVar = AccVar;

    % Initial State Vector
    X = [initialPosition';initialVelocity';initialOrientation'];

    % Initial State Covariance
    P = 0.1*eye(10);

    for i = 1 : length(accel)
        % IMU Readings
        w = gyro(i,:)';
        a = accel(i,:)';
        m = mag(i,:)';

        % State Prediction
        omega = Omega(w);
        q = X(7:end);
        dvdq = dxdq(a,q);
        F = eye(10) + [zeros(3)  ,eye(3)    ,zeros(3,4);
                       zeros(3)  ,zeros(3)  ,dvdq      ;
                       zeros(4,3),zeros(4,3),0.5*omega ] * delta;
        Rq = quat2rotm(q');
        Xp = [X(1:3) + X(4:6) * delta;
              X(4:6) + (Rq * a - [0,0,g]') * delta;
              X(7:end) + 0.5 * omega * delta * X(7:end)];

        % Covariance prediction
        W = 0.5 * delta * xi(q);
        posCov  = AccVar * delta^4 * eye(3);
        velCov  = AccVar * delta^2 * eye(3);
        RateCov = RateVar * (W) * (W');
        Q = [posCov    ,zeros(3)  ,zeros(3,4);
             zeros(3)  ,velCov    ,zeros(3,4)
             zeros(4,3),zeros(4,3),RateCov];     
        Pp = F * P * F' + Q;    

        % Measurement Matrix and Covariance
        qp = Xp(7:end);
        qpc = quatconj(qp')';
        gMeas = normalise(a);
        gMeas = globalG;
        mMeas = normalise(m);
        dgdq = dxdq(globalG,qpc);    
        dmdq = dxdq(globalM,qpc);
        C = quat2rotm(qpc');    
        GPSReady = ~mod(i-1, rateRatio);    
        z = [GPSPosition(floor((i-1)/rateRatio) + 1, :)';...
             GPSVelocity(floor((i-1)/rateRatio) + 1, :)';...
             gMeas;mMeas];
        h = [Xp(1:6); C*globalG; C*globalM];
        H = [eye(3)  ,zeros(3),zeros(3,4);
             zeros(3),eye(3)  ,zeros(3,4);
             zeros(3),zeros(3),dgdq      ;
             zeros(3),zeros(3),dmdq      ];
        R = [GPSPosVar*eye(3),zeros(3)        ,zeros(3)      ,zeros(3)     ;
             zeros(3)        ,GPSVelVar*eye(3),zeros(3)      ,zeros(3)     ;
             zeros(3)        ,zeros(3)        ,gravVar*eye(3),zeros(3)     ;
             zeros(3)        ,zeros(3)        ,zeros(3)      ,MagVar*eye(3)];    
       
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
        P = (eye(10) - K * H) * Pp;

        % Normalise Quaternion
        X(7:end) = normalise(X(7:end));

        % Save state vector.
        position(i, :) = X(1:3)';
        velocity(i, :) = X(4:6)';
        orientation(i, :) = X(7:end)';
    end
end