%% quaternionAttitudeEstimator.m
% Implements the initial attitude estimator described by section 7.5 in the
% report. This comes from the paper "Fast Quaternion Attitude Estimation
% from Two Vector Measurements" by Landis Markley
function q = quaternionAttitudeEstimator(b1, b2, r1, r2, a1, a2)
    % Normalise vectors
    b1 = normalise(b1);
    b2 = normalise(b2);
    r1 = normalise(r1);
    r2 = normalise(r2);
    
    % Compute cross products
    b3 = (cross(b1,b2))/norm(cross(b1,b2));
    r3 = (cross(r1,r2))/norm(cross(r1,r2));
    
    % Compute alpha, beta and gamma
    alpha = (1+dot(b3,r3))*(dot(a1*b1,r1)+dot(a2*b2,r2))+...
        dot(cross(b3,r3),cross(a1*b1,r1)+cross(a2*b2,r2));
    beta = dot((b3+r3),cross(a1*b1,r1)+cross(a2*b2,r2));
    gamma = sqrt(alpha^2+beta^2);
    
    % Compute Quaternion
    if alpha > 0
        q = (1/(2*sqrt(gamma*(gamma+alpha)*(1+dot(b3,r3))))) * ...
               [(gamma+alpha)*(1+dot(b3,r3)); ...
                (gamma+alpha)*cross(b3,r3)+beta*(b3+r3)];        
    else
        q = (1/(2*sqrt(gamma*(gamma+alpha)*(1+dot(b3,r3))))) * ...
               [beta*(1+dot(b3,r3)); ...
                beta*cross(b3,r3)+(gamma-alpha)*(b3+r3)];        
    end 
    
    % Reshape
    q = q';
end