% This example uses the Gauss-Newton algorithm to determine
% the minimum of Rosenbrock's loss function. This program 
% provides a plot of the function contours with the
% Gauss-Newton iterations.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example D.4

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Plot Space and Pre-Allocate Space
x1=[-3:0.05:3]';m1=length(x1);
x2=[-4:0.05:4]';m2=length(x2);
f=zeros(m1,m2);

% Get Function
ii=1;jj=1;
for i= -3:0.05:3,
  for j = -4:0.05:4,
    f(ii,jj)=100*(x2(jj)-x1(ii)^2)^2+(1-x1(ii))^2;
    jj=jj+1;
   end
ii=ii+1;
jj=1;
end

% Max Iterations, Initial Conditions
n_ite=4;
x0=[-1.2;1];xx=x0;
x=zeros(n_ite,2);x(1,:)=x0';

% Main Loop
for i = 1:n_ite,
 grad=[-400*x(i,1)*(x(i,2)-x(i,1)^2)-2*(1-x(i,1));200*(x(i,2)-x(i,1)^2)]; 
 s=-grad;
 hess=[-400*(x(i,2)-x(i,2)^2)+800*x(i,1)^2+2 -400*x(i,1);-400*x(i,1) 200];
 x(i+1,:)=x(i,:)+(inv(hess)*s)';
end
 
% Plot Results
clf
hold on
ccc=[2;10;50;100;150;200;250;300;350;400];
for i = 1:10,
 contour(x1,x2,f',[ccc(i) ccc(i)])
  if i == 1;
  axis([-3 3 -4 4])
 end
end
plot(x(:,1),x(:,2),'*')
plot(x(:,1),x(:,2))
hold off
set(gca,'fontsize',12)
xlabel('{\it x_1}')
ylabel('{\it x_2}')