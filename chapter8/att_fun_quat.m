function f = att_fun_quat(x,inertia,inertia_inv,u)

% Written by John L. Crassidis 9/03

% Function Routine for Attitude Dynamics

% Initialize
f=zeros(7,1);
q=x(1:4);q=q(:);
w=x(5:7);w=w(:);

% Matrices
qc=[0 -q(3) q(2);q(3) 0 -q(1);-q(2) q(1) 0];
pkiq=[q(4)*eye(3)+qc;-q(1:3)'];
wc=[0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0];

% Functions
f(1:4)=1/2*pkiq*w;
f(5:7)=-inertia_inv*wc*inertia*w+inertia_inv*u(:);