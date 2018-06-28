% This example shows a comparison between the extended Kalman
% filter and the Unscented filter to estimate the altitude,
% velocity and ballistic coefficient of a vertically falling
% body. This program provides a plot of the average magnitude
% of the position error by each filter using a Monte Carlo
% simulation consisting of 100 runs.

% Optimal Estimation of Dynamic Systems by (2nd ed.) Crassidis and Junkins
% Example 3.7

% Written by John L. Crassidis 9/03

% Other Required Routines: athansfun.m

% Initialize
dt=1;tf=60;t=[0:dt:tf]';m=length(t); % stores samples every second
runs=100;
xes1=zeros(m,runs);xes2=zeros(m,runs);xes3=zeros(m,runs);
xes1_kf=zeros(m,runs);xes2_kf=zeros(m,runs);xes3_kf=zeros(m,runs);

% Main Loop for Trial Runs
for jj = 1:runs,
    
    % Allocate Variables
    gam=5e-5;h=100000;mm=100000;
    x=zeros(m,3);x(1,1)=300000;x(1,2)=20000;x(1,3)=1e-3;
    ym=zeros(m,1);r=1e4;ym(1)=sqrt(mm^2+(x(1,1)-h)^2)+sqrt(r)*randn(1);
    
    % Initial Conditions
    pcov=diag([1e6 4e6 1e-4]);p=zeros(m,3);p(1,:)=diag(pcov)';
    xe=zeros(m,3);xe(1,1)=300000;xe(1,2)=20000;xe(1,3)=3e-5;
    
    % Kalman Filter Initial Conditions
    pcov_kf=diag([1e6 4e6 1e-4]);p_kf=zeros(m,3);p_kf(1,:)=diag(pcov_kf)';
    xe_kf=zeros(m,3);xe_kf(1,1)=300000;xe_kf(1,2)=20000;xe_kf(1,3)=3e-5;
    
    % Unscented Filter Parameters
    alp=1;beta=2;kap=0;n=3;
    lam=alp^2*(n+kap)-n;
    w0m=lam/(n+lam);
    w0c=lam/(n+lam)+(1-alp^2+beta);
    wim=1/(2*(n+lam));
    yez=zeros(1,6);
    
    % Interval Set to 1/64 Seconds for Integration
    dt=1/64;tf=60;t=[0:dt:tf]';m=length(t);
    % Main Filter Loop
    for i=1:m-1
        % Truth
        f1=dt*athansfun(x(i,:)',gam);
        f2=dt*athansfun(x(i,:)'+0.5*f1,gam);
        f3=dt*athansfun(x(i,:)'+0.5*f2,gam);
        f4=dt*athansfun(x(i,:)'+f3,gam);
        x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
        % Measurement
        ym(i+1)=sqrt(mm^2+(x(i+1,1)-h)^2)+sqrt(r)*randn(1);
    end
    x = x(1:1/dt:end,:);
    ym = ym(1:1/dt:end);
    
    % Interval Set to 1 Second for Estimate Propagation and Update
    dt=1;tf=60;t=[0:dt:tf]';m=length(t);
    for i=1:m-1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Unscented Filter
        % Covariance Decomposition
        psquare=chol(pcov)';
        sigv=real([sqrt(n+lam)*psquare -sqrt(n+lam)*psquare]);
        xx0=xe(i,:)';
        xx=sigv+kron(xe(i,:)',ones(1,2*n));
        
        % Calculate Mean Through Propagation
        f1=dt*athansfun([xx0 xx],gam);
        f2=dt*athansfun([xx0 xx]+0.5*f1,gam);
        f3=dt*athansfun([xx0 xx]+0.5*f2,gam);
        f4=dt*athansfun([xx0 xx]+f3,gam);
        xx0=xx0+1/6*(f1(:,1)+2*f2(:,1)+2*f3(:,1)+f4(:,1));
        xx=xx+1/6*(f1(:,2:2*n+1)+2*f2(:,2:2*n+1)+2*f3(:,2:2*n+1)+f4(:,2:2*n+1));
        
        xe(i+1,:)=w0m*xx0'+wim*sum(xx,2)';
        
        % Covariance
        pp0=w0c*(xx0-xe(i+1,:)')*(xx0-xe(i+1,:)')';
        pmat=xx-kron(xe(i+1,:)',ones(1,2*n));
        pcov=pp0+wim*pmat*pmat';
        
        % Output
        for j = 1:2*n
            yez(j)=sqrt(mm^2+(xx(1,j)-h)^2);
        end
        ye0=sqrt(mm^2+(xx0(1)-h)^2);
        
        ye=w0m*ye0+wim*sum(yez,2);
        
        % Calculate pyy
        pyy0=w0c*(ye0-ye)*(ye0-ye)';
        pyymat=yez-ye;
        pyy=pyy0+wim*pyymat*pyymat';
        
        % Calculate pxy
        pxy0=w0c*(xx0-xe(i+1,:)')*(ye0-ye);
        pxy=pxy0+wim*pmat*pyymat';
        
        % Innovations Covarinace
        pvv=pyy+r;
        
        % Gain and Update
        gain=real(pxy*inv(pvv));
        pcov=pcov-gain*pvv*gain';
        p(i+1,:)=diag(pcov)';
        xe(i+1,:)=xe(i+1,:)+(gain*(ym(i+1)-ye))';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Kalman Filter - Sometimes This Diverges
        % State Propagation
        f1=dt*athansfun(xe_kf(i,:)',gam);
        f2=dt*athansfun(xe_kf(i,:)'+0.5*f1,gam);
        f3=dt*athansfun(xe_kf(i,:)'+0.5*f2,gam);
        f4=dt*athansfun(xe_kf(i,:)'+f3,gam);
        xe_kf(i+1,:)=xe_kf(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
        
        % Estimate Output
        ye_kf=sqrt(mm^2+(xe_kf(i+1,1)-h)^2);
        
        % Covariance Propagation
        f=exp(-gam*xe_kf(i,1))*[0 -exp(gam*xe_kf(i,1)) 0
            gam*xe_kf(i,2)^2*xe_kf(i,3) -2*xe_kf(i,2)*xe_kf(i,3) -xe_kf(i,2)^2
            0 0 0];
        phi=c2d(f,zeros(3,1),dt);
        pcov_kf=phi*pcov_kf*phi';
        
        % Update
        h_kf=[(xe_kf(i+1,1)-h)/sqrt(mm^2+(xe_kf(i+1,1)-h)^2) 0 0];
        gain_kf=pcov_kf*h_kf'*inv(h_kf*pcov_kf*h_kf'+r);
        pcov_kf=(eye(3)-gain_kf*h_kf)*pcov_kf;
        p_kf(i+1,:)=diag(pcov_kf)';
        xe_kf(i+1,:)=xe_kf(i+1,:)+(gain_kf*(ym(i+1)-ye_kf))';
        
    end
    
    % Error for Each Trial Run
    xes1(:,jj)=abs(xe(:,1)-x(:,1));
    xes2(:,jj)=abs(xe(:,2)-x(:,2));
    xes3(:,jj)=abs(xe(:,3)-x(:,3));
    
    xes1_kf(:,jj)=abs(xe_kf(:,1)-x(:,1));
    xes2_kf(:,jj)=abs(xe_kf(:,2)-x(:,2));
    xes3_kf(:,jj)=abs(xe_kf(:,3)-x(:,3));
    
end

% Average and 3-Sigma Outlier
xerr=[sum(xes1,2) sum(xes2,2) sum(xes3,2)]/runs;
sig3=p.^(0.5)*3;
xerr_kf=[sum(xes1_kf,2) sum(xes2_kf,2) sum(xes3_kf,2)]/runs;
sig3_kf=p_kf.^(0.5)*3;

% Plot Results
plot(t,xerr(:,1),t,xerr_kf(:,1),'--')
set(gca,'Fontsize',12)
ylabel('Absolute Error of Average Altitude Error (M)')
xlabel('Time (Sec)')
legend('Unscented Filter','Extended Kalman Filter')