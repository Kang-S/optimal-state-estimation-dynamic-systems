% This example shows the power of weighted least squares 
% for a simple system. It outputs the current weight, the 
% estimated value and the constraint residual norm.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 1.3

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Measurements with Noise
t=[0:6/90:2]';m=length(t);
y=t+sin(t)+2*cos(2*t)-0.4*exp(t)/10000;
ym=y+0.1*randn(m,1);

% Set First Three Measurements to Truth
ym(1:3)=y(1:3);

% Basis Functions
h=[t sin(t) cos(2*t)];

% Weight Vector
weight=[1 1e1 1e2 1e5 1e7 1e10 1e15];

% Show Solutions with Different Weights
for i=1:length(weight)
  w=diag([weight(i) weight(i) weight(i) ones(1,28)]);
  current_w=weight(i)
  xe=(inv(h'*w*h)*h'*w*ym)'
  ye=xe(1)*t+xe(2)*sin(t)+xe(3)*cos(2*t);
  norm_error=norm(ye(1:3)-ym(1:3))
  disp(' Press any key to continue')
  disp(' ')
  pause
end