function [v,x_bound,y_bound]=ellipse_bound2(r,mean,m)
%function [v,x_bound,y_bound]=ellipse_bound2(r,mean,m)
%
% This m-file produces an m x 2 matrix of zero-mean 
% Gaussian noise with correlated covariance r and nonzero 
% mean. It also provides 3-sigma ellipse bounds.
%
%  The inputs are:
%       r = covariance matrix (2x2)
%    mean = mean (1x2 or 2x1)
%       m = number of points to produce
%
%  The outputs are:
%       v = noise matrix (mxn)
% x_bound = 3-sigma bound for x
% y_bound = 3-sigma bound for y

% Written by John L. Crassidis 1/4/06

% Get Eigenvalues and Angle
z=sqrt(r(1,1)^2+r(2,2)^2-2*r(1,1)*r(2,2)+4*r(1,2)^2);
lam1=(r(1,1)+r(2,2)-z)/2;
lam2=(r(1,1)+r(2,2)+z)/2;

% Get Angle of Rotation from Eigenvectors
den1=sqrt(r(1,2)^2+(lam1-r(1,1))^2);
den2=sqrt(r(1,2)^2+(lam2-r(1,1))^2);
if r(1,2) == 0 
 theta=0;
else
 theta=atan2(abs(r(1,2))/den2,r(1,2)/den1);
end

% Get Uncorrelated Noise
v_uncorr=randn(m,2)*diag([sqrt(lam1) sqrt(lam2)]); 

% Get Eigenvectors
u1=[r(1,2);lam1-r(1,1)]/den1;
u2=[r(1,2);lam2-r(1,1)]/den2;
u=[u1 u2];

% Get Correlated Noise
v=(u*v_uncorr')';
v(:,1)=v(:,1)+mean(1);
v(:,2)=v(:,2)+mean(2);

% Determine x and y Values
x3=sqrt(lam1)*3;
x_pos=[0:x3/50:x3]';
y_pos=3*(lam2*(1-x_pos.^2/x3^2)).^(0.5);
x=[x_pos;flipud(x_pos);-x_pos;-flipud(x_pos)]';
y=[y_pos;-flipud(y_pos);-y_pos;flipud(y_pos)]';

% Get Bounds Through Rotation
x_bound=x*cos(theta)+y*sin(theta)+mean(1);
y_bound=-x*sin(theta)+y*cos(theta)+mean(2);

% Plot Results
plot(x_bound,y_bound,v(:,1),v(:,2),'.')
