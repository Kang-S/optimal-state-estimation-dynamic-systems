% This example uses nonlinear least squares to determine
% the position of a vehicle on the Earth from GPS pseudorange
% measurements. It outputs the iterations and final 3-sigma 
% values.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 6.2

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% ECEF Position of DC
pos_dc=[1.132049606404574e+006;-4.903445558634290e+006;3.905453984115508e+006];

% SV Locations in ECEF Coordinates
sv1=[1.576473324923273e+007;-1.592675840044303e+006;2.124465542865848e+007];
sv2=[6.057534411078112e+006;-1.718695821956966e+007;1.939668903692545e+007];
sv3=[4.436748442200000e+006;-2.577117443480821e+007;-1.546041420776250e+006];
sv4=[-9.701586913814833e+006;-1.968746785735240e+007;1.535911814661023e+007];
sv5=[2.361749683191935e+007;-1.189936945888706e+007;1.492340517958122e+006];
sv6=[1.454007052476596e+007;-1.220196503121028e+007;1.835263217079338e+007];

% True Values and Starting Values
xt=[pos_dc;85000];
x=[0;0;0;0];

% Number of Iterations and Stored Values
nit=5;xe=zeros(nit+1,4);xe(1,:)=x';
sigm=5;

% Pseudorange Measurements
rho1=sqrt((sv1(1)-xt(1))^2 + (sv1(2)-xt(2))^2 + (sv1(3)-xt(3))^2) + xt(4)+sigm*randn(1);
rho2=sqrt((sv2(1)-xt(1))^2 + (sv2(2)-xt(2))^2 + (sv2(3)-xt(3))^2) + xt(4)+sigm*randn(1);
rho3=sqrt((sv3(1)-xt(1))^2 + (sv3(2)-xt(2))^2 + (sv3(3)-xt(3))^2) + xt(4)+sigm*randn(1);
rho4=sqrt((sv4(1)-xt(1))^2 + (sv4(2)-xt(2))^2 + (sv4(3)-xt(3))^2) + xt(4)+sigm*randn(1);
rho5=sqrt((sv5(1)-xt(1))^2 + (sv5(2)-xt(2))^2 + (sv5(3)-xt(3))^2) + xt(4)+sigm*randn(1);
rho6=sqrt((sv6(1)-xt(1))^2 + (sv6(2)-xt(2))^2 + (sv6(3)-xt(3))^2) + xt(4)+sigm*randn(1);

% Nonlinear Least Squares
for i = 1:nit,

 rho1e=sqrt((sv1(1)-x(1))^2 + (sv1(2)-x(2))^2 + (sv1(3)-x(3))^2) + x(4);
 rho2e=sqrt((sv2(1)-x(1))^2 + (sv2(2)-x(2))^2 + (sv2(3)-x(3))^2) + x(4);
 rho3e=sqrt((sv3(1)-x(1))^2 + (sv3(2)-x(2))^2 + (sv3(3)-x(3))^2) + x(4);
 rho4e=sqrt((sv4(1)-x(1))^2 + (sv4(2)-x(2))^2 + (sv4(3)-x(3))^2) + x(4);
 rho5e=sqrt((sv5(1)-x(1))^2 + (sv5(2)-x(2))^2 + (sv5(3)-x(3))^2) + x(4);
 rho6e=sqrt((sv6(1)-x(1))^2 + (sv6(2)-x(2))^2 + (sv6(3)-x(3))^2) + x(4);

 h1=[ -[(sv1(1)-x(1)) (sv1(2)-x(2)) (sv1(3)-x(3))]/(rho1e-x(4)) 1];
 h2=[ -[(sv2(1)-x(1)) (sv2(2)-x(2)) (sv2(3)-x(3))]/(rho2e-x(4)) 1];
 h3=[ -[(sv3(1)-x(1)) (sv3(2)-x(2)) (sv3(3)-x(3))]/(rho3e-x(4)) 1];
 h4=[ -[(sv4(1)-x(1)) (sv4(2)-x(2)) (sv4(3)-x(3))]/(rho4e-x(4)) 1];
 h5=[ -[(sv5(1)-x(1)) (sv5(2)-x(2)) (sv5(3)-x(3))]/(rho5e-x(4)) 1];
 h6=[ -[(sv6(1)-x(1)) (sv6(2)-x(2)) (sv6(3)-x(3))]/(rho6e-x(4)) 1];
 
 h=[h1;h2;h3;h4;h5;h6];

 dy=[rho1-rho1e;rho2-rho2e;rho3-rho3e;rho4-rho4e;rho5-rho5e;rho6-rho6e];
 dx=inv(h'*inv(sigm^2)*h)*h'*inv(sigm^2)*dy;

 x=x+dx;

 xe(i+1,:)=x';

 sig3=diag(inv(h'*inv(sigm^2)*h))'.^(0.5)*3;

end

% Show Results
iteration_results=xe
disp(' ')
sig3_outlier=sig3