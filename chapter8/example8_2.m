% This example show an LQR approach to asymptotically
% control a simple second-order system. This program 
% provides plots of the state trajectories.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 8.2

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% State Matrices, Time and Initial Conditions
a=[0 1;-2 2];b=[0;1];c=eye(2);d=zeros(2,1);
t=[0:0.01:10]';m=length(t);
x0=[5;3]';

% Weight Matrices and Gain Calculation
q=eye(2);r=0.1;
s=are(a,b*inv(r)*b',q);
gain=inv(r)*b'*s;

% Closed-Loop Dynamics Simulation
acl=a-b*gain;
x=lsim(acl,b,c,d,zeros(m,1),t,x0);

% Plot Results
subplot(211)
plot(t,x(:,1));
set(gca,'Fontsize',12);
axis([0 10 -2 6]);
set(gca,'Ytick',[-2 0 2 4 6])
xlabel('Time (Sec)')
ylabel('{\it x}_1({\it t})')

subplot(212)
plot(t,x(:,2));
set(gca,'Fontsize',12);
ylabel('{\it x}_2({\it t})')
xlabel('Time (Sec)')