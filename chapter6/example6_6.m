% This example uses ERA to identity the mass, stiffness and 
% damping matrices of a 4 mode system from simulated
% mass-position measurements. This program provides plots 
% of the simulated measurements and outputs the percent 
% damping, the mode singular values and MACs.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 6.6

% Written by John L. Crassidis 9/03

% Other Required Routines: era.m, hankl.m, p2m.m, svpm.m, deg2hz.m

% True Matrices
mm=eye(4);
kk=[10 -5 0 0;-5 10 -5 0;0 -5 10 -5; 0 0 -5 10];
cc=kk/5;

% State Space Model
a=[zeros(4) eye(4);-inv(mm)*kk -inv(mm)*cc];
b=[zeros(4);inv(mm)];
c=[eye(4) zeros(4)];
d=zeros(4,4);

% Truth
dt=0.1;
t=[0:dt:50]';
m=length(t);
u=[1;zeros(m-1,1)];

y1=lsim(a,b(:,1),c,d(:,1),u,t);
y2=lsim(a,b(:,2),c,d(:,2),u,t);
y3=lsim(a,b(:,3),c,d(:,3),u,t);
y4=lsim(a,b(:,4),c,d(:,4),u,t);

y=[y1 y2 y3 y4];

% Measurements
fac=0.01;
clear ym
sigm1=max(y(:,1))*fac;ym(:,1)=y(:,1)+sigm1*randn(m,1);
sigm2=max(y(:,2))*fac;ym(:,2)=y(:,2)+sigm2*randn(m,1);
sigm3=max(y(:,3))*fac;ym(:,3)=y(:,3)+sigm3*randn(m,1);
sigm4=max(y(:,4))*fac;ym(:,4)=y(:,4)+sigm4*randn(m,1);
sigm5=max(y(:,5))*fac;ym(:,5)=y(:,5)+sigm5*randn(m,1);
sigm6=max(y(:,6))*fac;ym(:,6)=y(:,6)+sigm6*randn(m,1);
sigm7=max(y(:,7))*fac;ym(:,7)=y(:,7)+sigm7*randn(m,1);
sigm8=max(y(:,8))*fac;ym(:,8)=y(:,8)+sigm8*randn(m,1);
sigm9=max(y(:,9))*fac;ym(:,9)=y(:,9)+sigm9*randn(m,1);
sigm10=max(y(:,10))*fac;ym(:,10)=y(:,10)+sigm10*randn(m,1);
sigm11=max(y(:,11))*fac;ym(:,11)=y(:,11)+sigm11*randn(m,1);
sigm12=max(y(:,12))*fac;ym(:,12)=y(:,12)+sigm12*randn(m,1);
sigm13=max(y(:,13))*fac;ym(:,13)=y(:,13)+sigm13*randn(m,1);
sigm14=max(y(:,14))*fac;ym(:,14)=y(:,14)+sigm14*randn(m,1);
sigm15=max(y(:,15))*fac;ym(:,15)=y(:,15)+sigm15*randn(m,1);
sigm16=max(y(:,16))*fac;ym(:,16)=y(:,16)+sigm16*randn(m,1);

% Call ERA
% Note: These Programs are Written By Juang and Pappa
[ad,bd,cd,dd,xs,xt,mac]=era(ym,4,4,8,100,0.1);

% Continuous-Time Results
[ae,be]=d2c(ad,bd,dt);
ce=cd;de=dd;
mass=inv(ce*ae*be);
ggg=-mass*ce*ae*ae*inv([ce;ce*ae]);
stiff=ggg(1:4,1:4);
damp=ggg(1:4,5:8);

% Plot Measurements
subplot(221)
plot(t,ym(:,1))
set(gca,'Fontsize',12);
grid
axis([0 50 -0.01 0.03]);
set(gca,'Xtick',[0 10 20 30 40 50]);
set(gca,'Ytick',[-0.01 0 0.01 0.02 0.03]);
xlabel('Time (Sec)')
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {\it {y}_1(t)} Measurement')

subplot(222)
plot(t,ym(:,2))
set(gca,'Fontsize',12);
grid
axis([0 50 -0.01 0.02]);
set(gca,'Xtick',[0 10 20 30 40 50]);
set(gca,'Ytick',[-0.01 0 0.01 0.02]);
xlabel('Time (Sec)')
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {\it {y}_2(t)} Measurement')

subplot(223)
plot(t,ym(:,3))
set(gca,'Fontsize',12);
grid
axis([0 50 -0.01 0.02]);
set(gca,'Xtick',[0 10 20 30 40 50]);
set(gca,'Ytick',[-0.01 0 0.01 0.02]);
xlabel('Time (Sec)')
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {\it {y}_3(t)} Measurement')

subplot(224)
plot(t,ym(:,1))
set(gca,'Fontsize',12);
grid
axis([0 50 -0.01 0.01]);
set(gca,'Xtick',[0 10 20 30 40 50]);
set(gca,'Ytick',[-0.01 -0.005 0 0.005 0.01]);
xlabel('Time (Sec)')
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {\it {y}_4(t)} Measurement')