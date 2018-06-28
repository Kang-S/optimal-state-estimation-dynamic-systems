% This example uses ERA to identity the mass, stiffness and 
% damping matrices of a 4 mode system from simulated
% mass-position measurements. The measurements are first 
% smoothed by an RTS smoother, which are inputted into the 
% ERA. This program provides plots of the position errors 
% with 3-sigma outliers. It also outputs the percent damping,
% the mode singular values and MACs.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 7.6

% Written by John L. Crassidis 9/03

% Other Required Routines: era.m, hankl.m, p2m.m, svpm.m, deg2hz.m

% True Matrices
mm=eye(4);
kk=0.9*[10 -5 0 0;-5 10 -5 0;0 -5 10 -5; 0 0 -5 10];
cc=kk/10;

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
fac=0.05;
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

% Process-Noise and Coversion to Discrete-Time
qcont=0.000001*eye(4);
gmat=[zeros(4);eye(4)];
bmat_ex=expm([-a gmat*qcont*gmat';zeros(8) a']*dt);
phi_dis=bmat_ex(9:16,9:16)';
qdis=phi_dis*bmat_ex(1:8,9:16);

% Discrete-Time State Matrices
[phi_dis,gam_dis]=c2d(a,b,dt);

% Measurement Covariances
rcov1=diag([sigm1^2 sigm2^2 sigm3^2 sigm4^2]);
rcov2=diag([sigm5^2 sigm6^2 sigm7^2 sigm8^2]);
rcov3=diag([sigm9^2 sigm10^2 sigm11^2 sigm12^2]);
rcov4=diag([sigm13^2 sigm14^2 sigm15^2 sigm16^2]);

% Covariances
pcov1=dare(phi_dis',c',qdis,rcov1,zeros(8,4),eye(8));
pcov2=dare(phi_dis',c',qdis,rcov2,zeros(8,4),eye(8));
pcov3=dare(phi_dis',c',qdis,rcov3,zeros(8,4),eye(8));
pcov4=dare(phi_dis',c',qdis,rcov4,zeros(8,4),eye(8));

% Gain Matrices
gain1=pcov1*c'*inv(c*pcov1*c'+rcov1);
gain2=pcov2*c'*inv(c*pcov2*c'+rcov2);
gain3=pcov3*c'*inv(c*pcov3*c'+rcov3);
gain4=pcov4*c'*inv(c*pcov4*c'+rcov4);

gain1b=(eye(8)-gain1*c)*pcov1*phi_dis'*inv(pcov1);
gain2b=(eye(8)-gain2*c)*pcov2*phi_dis'*inv(pcov2);
gain3b=(eye(8)-gain3*c)*pcov3*phi_dis'*inv(pcov3);
gain4b=(eye(8)-gain4*c)*pcov4*phi_dis'*inv(pcov4);

% Forward Kalman Filter
[ye1f,xe1f]=dlsim(phi_dis*(eye(8)-gain1*c),[gam_dis(:,1) phi_dis*gain1],c,[d(:,1) zeros(4,4)],[u ym(:,1:4)]);
[ye2f,xe2f]=dlsim(phi_dis*(eye(8)-gain2*c),[gam_dis(:,2) phi_dis*gain2],c,[d(:,2) zeros(4,4)],[u ym(:,5:8)]);
[ye3f,xe3f]=dlsim(phi_dis*(eye(8)-gain3*c),[gam_dis(:,3) phi_dis*gain3],c,[d(:,3) zeros(4,4)],[u ym(:,9:12)]);
[ye4f,xe4f]=dlsim(phi_dis*(eye(8)-gain4*c),[gam_dis(:,4) phi_dis*gain4],c,[d(:,4) zeros(4,4)],[u ym(:,13:16)]);

xe1fp=xe1f+(gain1*(ym(:,1:4)-xe1f(:,1:4))')';
xe2fp=xe2f+(gain2*(ym(:,5:8)-xe2f(:,1:4))')';
xe3fp=xe3f+(gain3*(ym(:,9:12)-xe3f(:,1:4))')';
xe4fp=xe4f+(gain4*(ym(:,13:16)-xe4f(:,1:4))')';

% Pre-Allocate Variables
xe1b=zeros(m,8);xe1b(m,:)=xe1f(m,:);
xe2b=zeros(m,8);xe2b(m,:)=xe2f(m,:);
xe3b=zeros(m,8);xe3b(m,:)=xe3f(m,:);
xe4b=zeros(m,8);xe4b(m,:)=xe4f(m,:);

% Backward Filter
for i = m-1:-1:1
 xe1b(i,:)=xe1fp(i,:)+(gain1b*(xe1b(i+1,:)'-xe1f(i+1,:)'))';
 xe2b(i,:)=xe2fp(i,:)+(gain2b*(xe2b(i+1,:)'-xe2f(i+1,:)'))';
 xe3b(i,:)=xe3fp(i,:)+(gain3b*(xe3b(i+1,:)'-xe3f(i+1,:)'))';
 xe4b(i,:)=xe4fp(i,:)+(gain4b*(xe4b(i+1,:)'-xe4f(i+1,:)'))';
end

ye1b=xe1b(:,1:4);
ye2b=xe2b(:,1:4);
ye3b=xe3b(:,1:4);
ye4b=xe4b(:,1:4);

% Smoother Covariance
pcov1s=dlyap(gain1b,(eye(8)-gain1*c)*pcov1-gain1b*pcov1*gain1b');
pcov2s=dlyap(gain2b,(eye(8)-gain2*c)*pcov2-gain2b*pcov2*gain2b');
pcov3s=dlyap(gain3b,(eye(8)-gain3*c)*pcov3-gain3b*pcov3*gain3b');
pcov4s=dlyap(gain4b,(eye(8)-gain4*c)*pcov4-gain4b*pcov4*gain4b');
sig1s=diag(pcov1s)'.^(0.5)*3;
sig2s=diag(pcov2s)'.^(0.5)*3;
sig3s=diag(pcov3s)'.^(0.5)*3;
sig4s=diag(pcov4s)'.^(0.5)*3;

% Call ERA
% Note: These Programs are Written By Juang and Pappa
[ad,bd,cd,dd,xs,xt,mac]=era([ye1b ye2b ye3b ye4b],4,4,8,100,0.1);

% Continuous-Time Results
[ae,be]=d2c(ad,bd,dt);
ce=cd;de=dd;
mass=inv(ce*ae*be);
ggg=-mass*ce*ae*ae*inv([ce;ce*ae]);
stiff=ggg(1:4,1:4);
damp=ggg(1:4,5:8);

% Plot Results
subplot(221)
plot(t,ye1b(:,1)-y(:,1),t,sig1s(1)*ones(m,1),t,-sig1s(1)*ones(m,1))
set(gca,'Fontsize',12);
grid
axis([0 50 -0.001 0.001]);
set(gca,'Xtick',[0 10 20 30 40 50]);
set(gca,'Ytick',[-0.001 -0.0005 0.0 0.0005 0.001]);
xlabel('Time (Sec)')
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {\it {y}_1(t)} Errors')

subplot(222)
plot(t,ye1b(:,2)-y(:,2),t,sig1s(2)*ones(m,1),t,-sig1s(2)*ones(m,1))
set(gca,'Fontsize',12);
grid
axis([0 50 -0.001 0.001]);
set(gca,'Xtick',[0 10 20 30 40 50]);
set(gca,'Ytick',[-0.001 -0.0005 0.0 0.0005 0.001]);
xlabel('Time (Sec)')
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {\it {y}_2(t)} Errors')

subplot(223)
plot(t,ye1b(:,3)-y(:,3),t,sig1s(3)*ones(m,1),t,-sig1s(3)*ones(m,1))
set(gca,'Fontsize',12);
grid
axis([0 50 -0.001 0.001]);
set(gca,'Xtick',[0 10 20 30 40 50]);
set(gca,'Ytick',[-0.001 -0.0005 0.0 0.0005 0.001]);
xlabel('Time (Sec)')
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {\it {y}_3(t)} Errors')

subplot(224)
plot(t,ye1b(:,4)-y(:,4),t,sig1s(4)*ones(m,1),t,-sig1s(4)*ones(m,1))
set(gca,'Fontsize',12);
grid
axis([0 50 -0.001 0.001]);
set(gca,'Xtick',[0 10 20 30 40 50]);
set(gca,'Ytick',[-0.001 -0.0005 0.0 0.0005 0.001]);
xlabel('Time (Sec)')
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {\it {y}_4(t)} Errors')