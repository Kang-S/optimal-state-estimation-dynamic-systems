% In this example the navigation Rao-Blackwellized Particle Filter (RBPF) 
% is used to track an unknown object's position and velocity using a set 
% of two range measurements. The states of the unknown object are its 
% planar position and associated velocity. 

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.14

% Written by John L. Crassidis 2/10

% Other Required Routines: correlated_noise, resample_pf

% Time
dt=0.1;tf=240*60;t=[0:dt:tf]';m=length(t);

% Positions
x_pos1=linspace(-5,30,m)';
y_pos1=zeros(m,1);

x_pos2=10*cos(0.001*t);
y_pos2=30*sin(0.005*t);

% Truth
phi=[1 0 dt 0;0 1 0 dt;0 0 1 0;0 0 0 1];
x0=[15 15 0 0]';x=zeros(m,4);x(1,:)=x0';
q=1e-10;
qd=q*[dt^3/3*eye(2) dt^2/2*eye(2);dt^2/2*eye(2) dt*eye(2)];
w=correlated_noise(qd,m);
for i =1:m-1, 
 x(i+1,:)=(phi*x(i,:)')'+w(i,:);
end

% Measurements
y=[((x_pos1-x(:,1)).^2+(y_pos1-x(:,2)).^2).^(0.5) ((x_pos2-x(:,1)).^2+(y_pos2-x(:,2)).^2).^(0.5)];
r=0.01*eye(2);
ym=y+sqrt(r(1,1))*randn(m,2);

% Particles
n=500;
w_particle=inv(n)*ones(n,1);
q1=qd(1:2,1:2);q2=qd(3:4,3:4);q12=qd(1:2,3:4);

% Initial Conditions
xe10=[0 0];xe20=[0 0];
x_particle1=kron(xe10,ones(n,1))+sqrt(64)*randn(n,2);
p0=diag([0.001 0.001]);
x_particle2=kron(xe20,ones(n,1))+sqrt(0.001)*randn(n,2);
p=p0;

xe=zeros(m,4);p_cov=zeros(m,4);

% Parameters for RBPF
c=q12'*pinv(q1);
d=eye(2)-dt*c;

% Main Loop
for i = 1:m
 res1=ym(i,1)-((x_pos1(i)-x_particle1(:,1)).^2+(y_pos1(i)-x_particle1(:,2)).^2).^(0.5);
 res2=ym(i,2)-((x_pos2(i)-x_particle1(:,1)).^2+(y_pos2(i)-x_particle1(:,2)).^2).^(0.5);
 w_particle=w_particle.*exp(-(res1.^2+res2.^2)/(2*r(1,1)));
 w_particle=w_particle/sum(w_particle);
 
 xe(i,:)=sum([w_particle.*x_particle1(:,1) w_particle.*x_particle1(:,2) w_particle.*x_particle2(:,1) w_particle.*x_particle2(:,2)]);
  
 % Resample
 [x_particle1,w_particle]=resample_pf(x_particle1,w_particle);
 f=x_particle1;
 
 noise_part1=correlated_noise(dt^2*p+q1,n);
 x_particle1=[f(:,1)+dt*x_particle2(:,1)+noise_part1(:,1) f(:,2)+dt*x_particle2(:,2)+noise_part1(:,2)];

 gain=dt*p*pinv(dt^2*p+q1);
 x_particle2=x_particle2+(gain*(x_particle1'-f'-dt*x_particle2'))';
 p=p-dt*gain*p;
 
 x_diff=[w_particle.*(x_particle1(:,1)-xe(i,1)) w_particle.*(x_particle1(:,2)-xe(i,2)) w_particle.*(x_particle2(:,1)-xe(i,3)) w_particle.*(x_particle2(:,2)-xe(i,4))];
 pcov=x_diff'*[x_particle1(:,1)-xe(i,1) x_particle1(:,2)-xe(i,2) x_particle2(:,1)-xe(i,3) x_particle2(:,2)-xe(i,4)]+[zeros(2) zeros(2);zeros(2) p];
 p_cov(i,:)=diag(pcov)';
 
 x_particle2=(d*x_particle2')'+(c*(x_particle1'-f'))';
 p=d*p*d'+q2-q12'*pinv(q1)*q12;
 
end  

% Plot Results
sig3=p_cov.^(0.5)*3;
plot(t/60,sig3(:,1),t/60,xe(:,1)-x(:,1),t/60,-sig3(:,1))
axis([0 240 -1 1])
set(gca,'XTick',[0 40 80 120 160 200 240])
set(gca,'fontsize',12)
xlabel('Time (Min)')
ylabel('Position Errors (km)')