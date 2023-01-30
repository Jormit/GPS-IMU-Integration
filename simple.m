%% simple.m
% Implements basic mechanization described by section 7.1.4.2 in the
% report.
function [velocity, position, orientation] = ... 
    simple(gyro, accel, ~, ~, ~, ...
    initialPosition, initialVelocity, initialOrientation, delta, ~, ~, ~)
    
    % Field Parameters
    g = [0,0,9.81];
    
    % Initialise      
    velocity = initialVelocity;
    position = initialPosition;
    orientation = initialOrientation;
    
    % Main algorithm
    for n = 1 : length(accel)     
        % Update orientation.
        w = [0,gyro(n,:)];
        w_dot = 0.5 * quatmultiply(orientation(n,:), w);
        orientation(n + 1,:) = normalise(orientation(n,:) + w_dot*delta);

        % Update position.
        a = rotatepoint(quaternion(orientation(n + 1,:)),accel(n,:)) - g;
        velocity(n + 1,:) = velocity(n,:) + a*delta;
        position(n + 1,:) = position(n,:) + velocity(n + 1,:)*delta;
    end
end