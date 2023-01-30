%% normalise.m
% Normalises a vector by dividing it by its norm.
function [xn] = normalise(x)
    xn = x / norm(x);
end
