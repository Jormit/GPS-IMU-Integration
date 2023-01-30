%% UKF.m
% Implements Uscented Kalman Filter described by section 8.4.4 and 7.3.5.10
% in the report.
function [velocity, position, orientation] = ...
    UKF(gyro, accel, mag, GPSPosition, GPSVelocity, ...
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

    % Initial State Vector
    X = [initialPosition';initialVelocity';initialOrientation'];

    % Initial State Covariance
    P = 0.1*eye(10);

    % Weighting Matrices
    alpha = 1e-3;
    ro = 0;
    beta = 2;
    n = length(X);
    lambda = alpha^2*(n + ro) - n;
    W_0_m = lambda/(n+lambda);
    W_0_c = lambda/(n+lambda) + (1 - alpha^2 + beta);
    W_i = 1/(2*(n+lambda));

    for i = 1 : length(accel)
        % IMU Readings
        w = gyro(i,:)';
        a = accel(i,:)';
        m = mag(i,:)';

        % State propagation
        X_sigma = zeros(n, 2*n + 1);
        X_sigma(:,1) = X;
        X_sigma(:,2:n+1) = X + chol(P*(n+lambda));
        X_sigma(:,n+2:2*n+1) = X - chol(P*(n+lambda));    
        Xp_sigma = f(X_sigma, w, a, [0,0,g]', delta);

        % State Prediction
        Xp = sum(W_i*Xp_sigma(:,2:2*n+1),2) + W_0_m*Xp_sigma(:,1);

        % Covariance Prediction
        q = X(7:end);
        W = 0.5 * delta * xi(q);
        posCov  = AccVar * delta^4 * eye(3);
        velCov  = AccVar * delta^2 * eye(3);
        RateCov = RateVar * (W) * (W');
        Q = [posCov    ,zeros(3)  ,zeros(3,4);
             zeros(3)  ,velCov    ,zeros(3,4)
             zeros(4,3),zeros(4,3),RateCov];  
        Pp = W_i*(Xp_sigma(:,2:2*n+1) - Xp)*(Xp_sigma(:,2:2*n+1) - Xp)' ... 
           + W_0_c*(Xp_sigma(:,1) - Xp)*(Xp_sigma(:,1) - Xp)' + Q;

        % Measurment Model propagation
        gMeas = normalise(a);
        gMeas = globalG;
        mMeas = normalise(m);    
        GPSReady = ~mod(i-1, rateRatio);
        Zp_sigma = h(Xp_sigma, globalG, globalM, GPSReady);
        if (GPSReady)
            Z = [GPSPosition((i-1)/rateRatio + 1, :)';...
                 GPSVelocity((i-1)/rateRatio + 1, :)';...
                 gMeas;mMeas];
            R = [GPSPosVar*eye(3),zeros(3)        ,zeros(3)     ,zeros(3)     ;
                 zeros(3)        ,GPSVelVar*eye(3),zeros(3)     ,zeros(3)     ;
                 zeros(3)        ,zeros(3)        ,AccVar*eye(3),zeros(3)     ;
                 zeros(3)        ,zeros(3)        ,zeros(3)     ,MagVar*eye(3)];        
        else
            Z = [gMeas; mMeas];
            R = [AccVar*eye(3),zeros(3)     ;
                 zeros(3)     ,MagVar*eye(3)];
        end

        % Measurment covariances
        Zp = W_i*sum(Zp_sigma(:,2:2*n+1),2) + W_0_m*Zp_sigma(:,1);
        Pzp = W_i*(Zp_sigma(:,2:2*n+1) - Zp)*(Zp_sigma(:,2:2*n+1) - Zp)'...
            + W_0_c*(Zp_sigma(:,1) - Zp)*(Zp_sigma(:,1) - Zp)';
        Pzpxp = W_i*(Xp_sigma(:,2:2*n+1) - Xp)*(Zp_sigma(:,2:2*n+1) - Zp)'...
              + W_0_c*(Xp_sigma(:,1) - Xp)*(Zp_sigma(:,1) - Zp)';

        % Update
        K = Pzpxp/(Pzp + R);
        X = Xp + K*(Z-Zp);
        P = Pp - K*Pzp*K';

        % Normalise Quaternion
        X(7:end) = normalise(X(7:end));    

        % Save state vector.
        position(i, :) = X(1:3)';
        velocity(i, :) = X(4:6)';
        orientation(i, :) = X(7:end)';    
    end
end

% Vectorised f and h
function Xp = f(X, rates, accel, gravity, delta)
    omega = Omega(rates);
    [~,cols] = size(X);
    Xp = zeros(size(X));
    for i = 1:cols        
        Rq = quat2rotm(X(7:end,i)');
        Xp(:,i) = [X(1:3,i) + X(4:6,i) * delta;
                   X(4:6,i) + (Rq * accel - gravity) * delta;
                   X(7:end,i) + 0.5 * omega * delta * X(7:end,i)];
    end
end

function Zp_sigma = h(Xp_sigma, globalG, globalM, GPSReady)
    [~,cols] = size(Xp_sigma);
    if (GPSReady)
        Zp_sigma = zeros(12,cols);
        for i = 1:cols   
            C = quat2rotm(Xp_sigma(7:end,i)')';
            Zp_sigma(:,i) = [Xp_sigma(1:6,i); C*globalG; C*globalM];
        end
    else
        Zp_sigma = zeros(6,cols);
        for i = 1:cols      
            C = quat2rotm(Xp_sigma(7:end,i)')';
            Zp_sigma(:,i) = [C*globalG; C*globalM];
        end
    end
end