% This example uses the parameters from example 1.2, but 
% solves the problem using nonlinear least squares. It 
% shows the iterations and final solutions converted 
% into discrete-time domain.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 1.7

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Get Truth and Measurements
dt=0.1;tf=10;
t=[0:dt:tf];
m=length(t);
a=-1;b=1;c=1;d=0;u=[100;zeros(m-1,1)];
[ad,bd]=c2d(a,b,dt);
y=dlsim(ad,bd,c,d,u);
ym=y+0.08*randn(m,1);
w=inv(0.08^2);

% Initialize Variables
h=[ym(1:m-1) u(1:m-1)];
clear xe
xe(1,:)=[5 5];
xx=100000;

% Nonlinear Least Squares
i=1;
while norm(xx) > 1e-8
if i > 50, break; end
aa=xe(i,1);bb=xe(i,2);
ea=exp(aa*dt);
h=[dt*ym(1:m-1)*ea+(bb/aa^2*(1-ea)+bb/aa*dt*ea)*u(1:m-1) 1/aa*(ea-1)*u(1:m-1)];
xx=inv(h'*w*h)*h'*w*(ym(2:m)-ym(1:m-1)*ea-bb/aa*(ea-1)*u(1:m-1));
xe(i+1,:)=xe(i,:)+xx';
i=i+1;
end

iteration_results=xe
disp(' ')
[phie,game]=c2d(xe(i,1),xe(i,2),dt)