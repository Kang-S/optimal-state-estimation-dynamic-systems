function out=dotprod1(x,y);
%function out=dotprod1(x,y);
%
% This program takes the dot product of
% each x(i,:) with each y(i,:). The
% output is a mx1 vector.

% Written by John L. Crassidis 4/24/95

col=size(x,2);
xs=x.*y;
out=cumsum(xs')';
out=out(:,col);
