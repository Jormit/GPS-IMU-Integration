%% Omega.m
% Generates the Omega matrix described by section 7.1.4.2 in the report.
function Omega = Omega(x)
    Omega = [0   , -x(1), -x(2), -x(3);
             x(1),  0   ,  x(3), -x(2);
             x(2), -x(3),  0   ,  x(1);
             x(3),  x(2), -x(1),  0; ];
end
