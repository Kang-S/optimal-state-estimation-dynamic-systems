% This plots Himmelblau's function as shown by Figure C.1.
% It also plots the gradient method iterations for four 
% initial guesses.

% Optimal Estimation of Dynamic Systems by Crassidis and Junkins
% Figure C.1

% Written by John L. Crassidis 9/03

% Other Required Routines: himmelblau_fun.m

% Set Parameter Space
x1=[-6:0.1:6]';m1=length(x1);
x2=[-6:0.1:6]';m2=length(x2);
f=zeros(m1,m2);

% Get Points from Himmelblau's Function
ii=1;jj=1;
for i= -6:0.1:6,
  for j = -6:0.1:6,
    f(ii,jj)=(x1(ii)^2+x2(jj)-11)^2+(x1(ii)+x2(jj)^2-7)^2;
    jj=jj+1;
   end
ii=ii+1;
jj=1;
end

% Number of Iterations of Gradient Method
n_ite=5;

% Initial Starting Points
x0_1=[-2;5];
x0_2=[-2;-5];
x0_3=[2;5];
x0_4=[5;-5];
xx=[x0_1';x0_2';x0_3';x0_4'];

% Loop for 4 Points
for j = 1:4,

 x=zeros(n_ite,2);x(1,:)=xx(j,:);
 alp0=3;alpp=alp0;
 alp=zeros(n_ite,1);alp(1)=alp0;    

% Gradient Method 
 for i = 1:n_ite,
  grad=[4*x(i,1)*(x(i,1)^2+x(i,2)-11)+2*(x(i,1)+x(i,2)^2-7)
       2*(x(i,1)^2+x(i,2)-11)+4*x(i,2)*(x(i,1)+x(i,2)^2-7)]; 
   s=-grad;
   alp(i)=fminsearch('himmelblau_fun',alp(i),[],s,x(i,:)');
   x(i+1,:)=x(i,:)+alp(i)*s';
 end

% Plot Contours and Iterations 
 if j == 1,
  clf
  hold on
  ccc=[5;20;35;50;65;80;95;110;125;140];
  for i = 1:10,
   contour(x1,x2,f',[ccc(i) ccc(i)])
   if i == 1;
    axis([-6 6 -6 6])
  end
  end
 end
 plot(x(:,1),x(:,2),'*')
 plot(x(:,1),x(:,2))
end
hold off

% Fontsize and Labels
set(gca,'fontsize',12)
xlabel('{\it x_1}')
ylabel('{\it x_2}')
%xlabel('x1')
%ylabel('x2')