% This example shows a simplified version the simultaneous localization 
% and mapping (SLAM) problem, where the platform motion and feature point 
% locations are two dimensional in nature.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 6.3

% Written by Manoranjan Majji 03/11

% Other Required Routines: none

% Sine Wave Profile for Triangulation
Npts = 50;
x = linspace(0, 10, Npts); 
y = 15 + 5*sin(x/2);

% Location 1
x1 = 1; y1 = 10;

pt1(1,:) = 1 - linspace(-0.5,0.5,10);
pt1(2,:) = 10*ones(1,10);

% Translation as Expressed in the WCS
t1hat(1,1) = -x1; t1hat(2,1) = -y1;

% Measurements
ytil(1,:) = x(1,5:3:20);
ytil(2,:) = y(1,5:3:20);

m = length(ytil(1,:));

% Location 2
x2 = 5; y2 = 15;

% Translation as Expressed in the WCS
t2hat(1,1) = -(x2 - x1); t2hat(2,1) = -(y2 - y1);

% Seen in the First Camera
ytil1 = ytil + t1hat*ones(1,m);

th2 = pi/4;
R2(1,1) = cos(th2); R2(1,2) = sin(th2); R2(2,1) = -R2(1,2);
R2(2,2) = cos(th2);

% Seen in the Second Camera
ytil2 = R2*(ytil1 + t2hat*ones(1,m)); 

pt2 = (R2'*(pt1 - [1 ; 10]*ones(1,10)) + [5 ; 15]*ones(1,10)) ;

dtil = ytil2 - ytil1;
stil = ytil2 + ytil1;

% Least Squres Solution
for k = 1:m    
    err(2*k-1:2*k,1) = dtil(:,k);
    
    Htil(2*k-1, 1) = stil(2,k);
    Htil(2*k-1, 2:3) = [1, 0];
    Htil(2*k, 1) = -stil(1,k);
    Htil(2*k, 2:3) = [0, 1];
end

xhat = pinv(Htil'*Htil)*Htil'*err

q3hat = xhat(1);
IPQ = [1 -q3hat; q3hat, 1];

Rhat = inv(IPQ)*[1 q3hat; -q3hat, 1]

that = inv(IPQ)*xhat(2:3,1)

t2hat1 = -Rhat'*that;

yerr = err - Htil*xhat

% Plot Results
plot(x, y, ytil(1,:), ytil(2,:), '*', pt1(1,:),pt1(2,:),'k.',pt2(1,:),pt2(2,:),'k.',x1, y1, 'ko',x2, y2, 'ko');
axis([0 10 10 20]);
set(gca,'fontsize',12)
xlabel('x'); ylabel('y');
