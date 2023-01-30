%% Omega.m
% Generates the Xi matrix described by section 7.1.4.2 in the report.
function [xi] = xi(q)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    xi = [-q1, -q2, -q3;...
           q0, -q3,  q2;...
           q3,  q0, -q1;...
          -q2,  q1,  q0];
end

