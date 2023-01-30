%% enu2geo.m
% Wrapper around enu2geodetic which converts between ENU and Geodetic 
% coordinate systems.
function [positionLLA] = enu2geo(position, referenceLocation, ellipsoid)
    [lat,long,h] = enu2geodetic(...
        position(:,1),position(:,2),position(:,3),...
        referenceLocation(1),referenceLocation(2),referenceLocation(3),ellipsoid);
    positionLLA = [lat,long,h];    
end