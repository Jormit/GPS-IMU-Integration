%% Omega.m
% Generates jacobian of a 3d vector with respect to a quaternion described
% by section 8.4.3 in the report.
function [dxdq] = dxdq(x,q)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    x1 = x(1); x2 = x(2); x3 = x(3);
    dxdq = 2*...
    [q2*x3-q3*x2,-q2*x2-q3*x3,2*q2*x1-q1*x2-q0*x3,q0*x2-q1*x3+2*q3*x1;...
     q3*x1-q1*x3,q0*x3+2*q1*x2-q2*x1,-q1*x1-q3*x3,2*q3*x2-q2*x3-q0*x1;...
     q1*x2-q2*x1,2*q1*x3-q0*x2-q3*x1,q0*x1+2*q2*x3-q3*x2,-q1*x1-q2*x2];
end