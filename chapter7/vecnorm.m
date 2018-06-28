function out=vecnorm(x);
%function out=vecnorm(x);
%
% This program takes the norm of each 
% x(i,:) to produce an vector norm output, i.e.
% out = [x(:,1).^2 + x(:,2).^2 + ... ].^0.5

% Written by John L. Crassidis 4/22/95

col=size(x,2);
xs=x.*x;
out=cumsum(xs')';
out=out(:,col).^0.5;