% This example uses the method of gradients to determine the
% minimum of a quadratic function. This program provides a 
% plot of the function contours with the gradient iterations.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example D.3

% Written by John L. Crassidis 9/03

% Other Required Routines: quadratic_fun.m

% Plot Space and Pre-Allocate Space
x1=[-3:0.05:3]';m1=length(x1);
x2=[-3:0.05:3]';m2=length(x2);
f=zeros(m1,m2);

% Get Function
ii=1;jj=1;
for i= -3:0.05:3,
  for j = -3:0.05:3,
    f(ii,jj)=4*x1(ii)^2+3*x2(jj)^2-4*x1(ii)*x2(jj)+x1(ii);
    jj=jj+1;
   end
ii=ii+1;
jj=1;
end

% Max Iterations, Initial Conditions and 1-D Search Parameter
n_ite=10;
x0=[-1;3];xx=x0;
x=zeros(n_ite,2);x(1,:)=x0';
alp0=3;alpp=alp0;
alp=zeros(n_ite,1);alp(1)=alp0;

% Main Loop
for i = 1:n_ite,
 grad=[8*x(i,1)-4*x(i,2)+1;6*x(i,2)-4*x(i,1)]; 
 s=-grad;
 alp(i)=fminsearch('quadratic_fun',alp(i),[],s,x(i,:)');
 x(i+1,:)=x(i,:)+alp(i)*s';
end

% Plot Results
clf
hold on
ccc=[0.5;5.5;10.5;15.5;20.5;25.5;30.5;35.5;40.5;45.5];
for i = 1:10,
 contour(x1,x2,f',[ccc(i) ccc(i)])
 if i == 1;
  axis([-3 3 -3 3])
 end
end
plot(x(:,1),x(:,2),'*')
plot(x(:,1),x(:,2))
hold off
set(gca,'fontsize',12)
xlabel('{\it x_1}')
ylabel('{\it x_2}')