% This examples provides a plot of the inefficiency for 
% varying parameters in a simple system.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 2.10

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Set R and Clear cost Variable
r=eye(2);clear cost
r=[1 0; 0 1];

% Get Mesh
jj=1;
ii=1;
for i=1:-0.1:-.99
 for j=1:-0.1:-.99
  rp=[1+i 0; 0 1+j];
  s=svd(rp^(-0.5)*r*rp^(-0.5));
  lmax=max(s);
  lmin=min(s);
  e=(lmax+lmin)^2/4/lmax/lmin;
  cost(ii,jj)=e;
  jj=jj+1;
 end
 ii=ii+1;
 jj=1;
end

% Plot Results
x=[1:-0.1:-.99]';
y=[1:-0.1:-.99]';
mesh(x,y,cost);
set(gca,'fontsize',12);
set(gca,'xtick',[-1 -0.5 0 0.5 1]);
set(gca,'ytick',[-1 -0.5 0 0.5 1]);
set(gca,'ztick',[0 1 2 3 4 5 6]);
xlabel('{\alpha}')
ylabel('{\beta}')
zlabel('Inefficiency Bound')