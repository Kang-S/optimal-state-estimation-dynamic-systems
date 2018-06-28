% This example shows a plot of the forward filter, backward 
% filter and smoother covariances for a simple continuous-time 
% first-order system.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 5.2

% Written by John L. Crassidis 9/03

% Other Required Routines: ric_forfun.m ric_backfun.m

% Time Vectors
t=[0:0.1:10];
tt=[10:-0.1:0]';

% Forwards and Backwards Integrations
[t,pf]=ode23(@ric_forfun,t,1000);
[t,pbi]=ode23(@ric_backfun,t,0);

% Covariances
pb=pbi.^(-1);pb(1)=1000;
p=(pf.^(-1)+pbi).^(-1);

% Plot R=esults
semilogy(tt,pb,'--',t,pf,'--',t,p)
set(gca,'Fontsize',12)
ylabel('Covariances')
xlabel('Time (Sec)')
set(gca,'xtick',[0 1 2 3 4 5 6 7 8 9 10])
text(3.8+0.4,3.5,'Backward Filter','Fontsize',12)
text(3.9+0.4,0.93,'Forward Filter','Fontsize',12)
text(4.1+0.4,0.47,'Smoother','Fontsize',12)