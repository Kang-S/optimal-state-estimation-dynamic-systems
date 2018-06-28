% This example uses the parameters from example 1.3, but 
% solves the problem using constrained least squares. 
% Three cases are shown: 1) first measurement is set to the 
% true value, 2) first two measurements are set to the true
% values, and 3) first three measurements are set to the 
% true values. Results for all three cases are shown.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 1.4

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Measurements with Noise
t=[0:6/90:2]';m=length(t);
y=t+sin(t)+2*cos(2*t)-0.4*exp(t)/10000;
ym=y+0.1*randn(m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 1
% Set First Measurement to Truth
ym(1:1)=y(1:1);
ym1=ym(2:31);
ym2=ym(1);

% Basis Functions
h1=[t(2:31) sin(t(2:31)) cos(2*t(2:31))];
h2=[t(1:1) sin(t(1:1)) cos(2*t(1:1))];

% Get Constrained Estimate
xe_bar_case1=(inv(h1'*h1)*h1'*ym1)'
gain=inv(h1'*h1)*h2'*inv(h2*inv(h1'*h1)*h2');
xe_case1=(xe_bar_case1'+gain*(ym2-h2*xe_bar_case1'))'

disp(' Press any key to continue')
disp(' ')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 2
% Set First Two Measurements to Truth
ym(1:2)=y(1:2);
ym1=ym(3:31);
ym2=ym(1:2);

% Basis Functions
h1=[t(3:31) sin(t(3:31)) cos(2*t(3:31))];
h2=[t(1:2) sin(t(1:2)) cos(2*t(1:2))];

% Get Constrained Estimate
xe_bar_case2=(inv(h1'*h1)*h1'*ym1)'
gain=inv(h1'*h1)*h2'*inv(h2*inv(h1'*h1)*h2');
xe_case2=(xe_bar_case2'+gain*(ym2-h2*xe_bar_case2'))'

disp(' Press any key to continue')
disp(' ')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 3
% Set First Three Measurements to Truth
ym(1:3)=y(1:3);
ym1=ym(4:31);
ym2=ym(1:3);

% Basis Functions
h1=[t(4:31) sin(t(4:31)) cos(2*t(4:31))];
h2=[t(1:3) sin(t(1:3)) cos(2*t(1:3))];

% Get Constrained Estimate
xe_bar_case3=(inv(h1'*h1)*h1'*ym1)'
gain=inv(h1'*h1)*h2'*inv(h2*inv(h1'*h1)*h2');
xe_case3=(xe_bar_case3'+gain*(ym2-h2*xe_bar_case3'))'