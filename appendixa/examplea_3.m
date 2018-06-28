% This example shows a linear perturbation technique to 
% study the behavior of a highly maneuverable aircraft 
% which exhibits nonlinear behavior. This program provides 
% plots of the state perturbation trajectories and final 
% state trajectories.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example A.3

% Written by John L. Crassidis 9/03

% Other Required Routines: f8_fun.m

% Initialize Variables
m=501;dt=.01;
t=[0:dt:m*dt-dt]';
x0=[25*pi/180;0;0];
x=zeros(m,3);x(1,:)=x0';
xn=zeros(m,3);xn(1,:)=[x0(1)-1*pi/180;0;0]';
dx=zeros(m,3);dx(1,:)=[1*pi/180;0;0]';

% True and Nominal Values
c=[1;1;0.09;0.88;0.47;3.85;0.01;.396;4.208;0.47;3.564];
cn=c*1;

% Main Loop for Integration
for i=1:m-1,

 f1=dt*f8_fun(x(i,:),c);
 f2=dt*f8_fun(x(i,:)+0.5*f1',c);
 f3=dt*f8_fun(x(i,:)+0.5*f2',c);
 f4=dt*f8_fun(x(i,:)+f3',c);
 x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');

 f1=dt*f8_fun(xn(i,:),cn);
 f2=dt*f8_fun(xn(i,:)+0.5*f1',cn);
 f3=dt*f8_fun(xn(i,:)+0.5*f2',cn);
 f4=dt*f8_fun(xn(i,:)+f3',cn);
 xn(i+1,:)=xn(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');

 x1=xn(i,1);x2=xn(i,2);x3=xn(i,3);
 a11=-2*c(2)*x1*x3-c(3)*x3-c(4)+2*c(5)*x1+3*c(6)*x1*x1;
 a12=-2*c(7)*x2;
 a13=c(1)-c(2)*x1^2-c(3)*x1;
 a21=0;a22=0;a23=1;
 a31=-c(9)-2*c(10)*x1-3*c(11)*x1^2;
 a32=0;
 a33=-c(8);
 a=[a11 a12 a13;a21 a22 a23;a31 a32 a33];
 
 phi=c2d(a,zeros(3,1),dt);
 dx(i+1,:)=(phi*dx(i,:)')';
 
end

%Plot Results
plot(t,dx(:,1)*180/pi,t,dx(:,2)*180/pi,'--',t,dx(:,3)*180/pi,'-.');
set(gca,'Fontsize',12);
ylabel('State Perturbations (Deg)')
xlabel('Time (Sec)');
ax=legend('{{\it \delta}{\it x}_1}','{{\it \delta}{\it x}_2}','{{\it \delta}{\it x}_3}');
leg=findobj(ax,'type','text');
set(leg,'FontUnits','points','fontsize',12);
grid

disp(' Press any key to continue')
pause

plot(t,x(:,1)*180/pi,t,x(:,2)*180/pi,'--',t,x(:,3)*180/pi,'-.');
set(gca,'Fontsize',12);
ylabel('States (Deg)')
xlabel('Time (Sec)');
ax=legend('{{\it x}_1}','{{\it x}_2}','{{\it x}_3}');
leg=findobj(ax,'type','text');
set(leg,'FontUnits','points','fontsize',12);
grid