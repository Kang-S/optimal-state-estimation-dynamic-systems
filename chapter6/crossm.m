function out=crossm(x)
%function out=crossm(x)
%
% This program forms the (3x3) cross matrix:
%               [0 -x(3) x(2)]
%         out = [x(3) 0 -x(1)]
%               [-x(2) x(1) 0]

% Written by John L. Crassidis 4/24/95

x=x(:);
out=[0 -x(3) x(2);x(3) 0 -x(1);-x(2) x(1) 0];