%% enu2geo.m
% Wrapper around geodetic2enu which converts between Geodetic and ENU 
% coordinate systems.
function [positionEnu] = geo2enu(position, referenceLocation, ellipsoid)
    [x,y,z] = geodetic2enu(position(:,1), position(:,2), position(:,3), ...
        referenceLocation(1), referenceLocation(2), referenceLocation(3), ellipsoid);
    positionEnu = [x,y,z];
end