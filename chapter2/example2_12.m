% This example uses the parameters from example 1.2, but solves 
% the problem using total least squares. It outputs the current 
% value for the sampling interval, and bias estimates and mean square 
% errors using both standard least squares and total least squares.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 2.12

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Loop for Different dt
dt_vec=[2;1;0.5;0.1;0.05;0.01;0.005;0.001];
for j=1:length(dt_vec);

% Get Truth
 dt=dt_vec(j);tf=10;
 t=[0:dt:tf];
 m=length(t);
 y=zeros(m,1);
 a=-1;b=1;c=1;d=0;u=[100*.1/dt;zeros(m-1,1)];%u=ones(m,1);
 [ad,bd]=c2d(a,b,dt);
 y=dlsim(ad,bd,c,d,u);

% Simulation Parameters
 mm=1000;
 xls=zeros(mm,2);
 xtls=zeros(mm,2);

% Main Loop
 for i=1:mm,
  ym=y+0.08*randn(m,1);

% Least Squares Solution 
  w=inv(0.08^2);
  h=[ym(1:m-1) u(1:m-1)];
  x=inv(h'*h)*h'*ym(2:m);
  xls(i,:)=x'; 

% Total Least Squares Solution
  [uu,ss,vv]=svd([h ym(2:m)],0);
  xx=-inv(vv(3,3))*vv(1:2,3);
  xtls(i,:)=xx';

 end

% Show Results
 lse_ad=xls(:,1)-ad;lse_ad=(lse_ad.*lse_ad).^(0.5);
 tlse_ad=xtls(:,1)-ad;tlse_ad=(tlse_ad.*tlse_ad).^(0.5);

 lse_bd=xls(:,2)-bd;lse_bd=(lse_bd.*lse_bd).^(0.5);
 tlse_bd=xtls(:,2)-bd;tlse_bd=(tlse_bd.*tlse_bd).^(0.5);

 abs(mean(xls)-[ad bd])./[ad bd]*100;
 abs(mean(xtls)-[ad bd])./[ad bd]*100;

 bias_ls=[mean(xls(:,1))-ad mean(xls(:,2))-bd]'...
        *[mean(xls(:,1))-ad mean(xls(:,2))-bd];
 bias_tls=[mean(xtls(:,1))-ad mean(xtls(:,2))-bd]'...
         *[mean(xtls(:,1))-ad mean(xtls(:,2))-bd];

 cov_ls=cov(xls);
 cov_tls=cov(xtls);
 
 dt_current=dt
 
 bias_ls2=diag(bias_ls)'.^(0.5)
 bias_tls2=diag(bias_tls)'.^(0.5)

 mse_ls2=diag(cov_ls+bias_ls)'.^(0.5)
 mse_tls2=diag(cov_tls+bias_tls)'.^(0.5)
 
 disp(' Press any key to continue')
 disp(' ')
 pause
end


