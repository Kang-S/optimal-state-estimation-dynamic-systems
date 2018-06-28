% This example shows the usefulness of the Loop Transfer Recovery 
% approach for a simple continuous-time second-order system.
% It outputs the original phase margin, the trial weight, the 
% steady-state covariance and the current phase margin.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 8.3

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% State Model and Original Weights
f=[0 1;-3 -4];b=[0;1];h=[2 1];g=[35;-61];
num_g=[1 2];den_g=[1 4 3];
r=1;q=1;

% Gain and Trial Weights
l_gain=[50 10];
q2=[0;100;500;1e3;1e4];
m=length(q2);

% Place Estimator Poles at -22+-17.86j
k_gain=(f^2+44*f+802.9796*eye(2))*inv([h;h*f])*[0;1];
[num_k,den_k]=ss2tf(f-b*l_gain-k_gain*h,k_gain,l_gain,0,1);
num=conv(num_k,num_g);den=conv(den_k,den_g);

% Original Phase Margin
[gm,ph]=margin(num,den);
phase_orginal=ph

disp(' Press any key to continue')
disp(' ')
pause

% Main Loop for All Weights
for i = 1:m

% Riccati Solution and Gain
 p=are(f',h'*inv(r)*h,g*q*g'+q2(i)*[0;1]*[0 1]);
 k_gain=p*h'*inv(r);
 
% Get New Phase Margin
 [num_k,den_k]=ss2tf(f-b*l_gain-k_gain*h,k_gain,l_gain,0,1);
 num=conv(num_k,num_g);den=conv(den_k,den_g);
 [gm,ph]=margin(num,den);
 
% Show Results
 q_square=q2(i)
 p_steady=p
 gain=k_gain
 current_phase = ph
 
 disp(' Press any key to continue')
 disp(' ')
 pause
 
end