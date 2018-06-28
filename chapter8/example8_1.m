% This example shows the control of a rigid body constrained 
% to rotate about a fixed axis. This program provides plots 
% of the optimal rest-to-rest maneuver, and the effect of 
% the final time variation for a spinup maneuver.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 8.1

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Time, Control Input, Angle and Rate
t=[0:0.01:1]';
u=3*pi*(1-2*t);
phi=3*pi*(t.*t/2-t.*t.*t/3);
phid=3*pi*(t-t.*t);

% Plot Results
subplot(311)
plot(t,phi*180/pi)
set(gca,'fontsize',12);
ylabel('Angle (Deg)')
axis([0 1 0 90])
set(gca,'xtick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
set(gca,'ytick',[0 30 60 90])

subplot(312)
plot(t,phid*180/pi)
set(gca,'fontsize',12);
ylabel('Velocity (Deg/Sec)')
axis([0 1 0 180])
set(gca,'xtick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
set(gca,'ytick',[0 60 120 180])

subplot(313)
plot(t,u)
set(gca,'fontsize',12);
ylabel('Control Input')
xlabel('Time (Sec)')
set(gca,'xtick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
set(gca,'ytick',[-10 -5 0 5 10])

disp(' Press any key to continue')
pause
clf

% Case 1
tf=3*pi/2;
theta_f=pi/2;theta_0=0;thetad_f=1;thetad_0=0;
a1=theta_0;
a2=thetad_0;
a3=3*(theta_f-theta_0)/tf^2-(2*thetad_0+thetad_f)/tf;
a4=-2*(theta_f-theta_0)/tf^3+(thetad_0+thetad_f)/tf^2;
t1=[0:0.01:tf]';
u1=2*a3+6*a4*t1;
theta1=a1+a2*t1+a3*t1.^2+a4*t1.^3;

% Case 2
tf=3*pi/2+1;
theta_f=pi/2;theta_0=0;thetad_f=1;thetad_0=0;
a1=theta_0;
a2=thetad_0;
a3=3*(theta_f-theta_0)/tf^2-(2*thetad_0+thetad_f)/tf;
a4=-2*(theta_f-theta_0)/tf^3+(thetad_0+thetad_f)/tf^2;
t2=[0:0.01:tf]';
u2=2*a3+6*a4*t2;
theta2=a1+a2*t2+a3*t2.^2+a4*t2.^3;

% Case 3
tf=3*pi/2-1;
theta_f=pi/2;theta_0=0;thetad_f=1;thetad_0=0;
a1=theta_0;
a2=thetad_0;
a3=3*(theta_f-theta_0)/tf^2-(2*thetad_0+thetad_f)/tf;
a4=-2*(theta_f-theta_0)/tf^3+(thetad_0+thetad_f)/tf^2;
t3=[0:0.01:tf]';
u3=2*a3+6*a4*t3;
theta3=a1+a2*t3+a3*t3.^2+a4*t3.^3;

% Plot Results
subplot(211)
plot(t1,u1,t2,u2,'--',t3,u3,'-.')
axis([0 6 -0.1 0.5])
set(gca,'fontsize',12);
set(gca,'xtick',[0 1 2 3 4 5 6])
set(gca,'ytick',[-0.1 0 0.1 0.2 0.3 0.4 0.5])
legend('Case 1', 'Case 2','Case 3',4)
ylabel('Control Input')
xlabel('Time (Sec)')

subplot(212)
plot(t1,theta1*180/pi,t2,theta2*180/pi,'--',t3,theta3*180/pi,'-.')
axis([0 6 -10 90])
set(gca,'fontsize',12);
set(gca,'xtick',[0 1 2 3 4 5 6])
set(gca,'ytick',[-10 10 30 50 70 90])
legend('Case 1', 'Case 2','Case 3',2)
ylabel('Angle (Deg)')
xlabel('Time (Sec)')
hold on
tt=[0:1:6]';
plot(tt,zeros(7,1),':')
hold off