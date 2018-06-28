function q=qmethod(im,bm,wi,av,flag);
%function q=qmethod(im,bm,wi,av,flag);
%
% This program computes attitudes using the Q method.
% It also tries to make the quaternions continuous.
% This must have at least two vectors measurements at any time.
%
%  The inputs are:
%     im = inertial measurements [mx(3*s)], s = number of sensors
%     bm = body measurements  [mx(3*s)]
%     wi = weighting vector (sx1)
%     av = 1 for sensor available, 0 for not  [mxs]
%     al = alignment matrices for each sensor [3x(3*s)]
%   flag = 1 to try to make the quaternions continuous

% Written by John L. Crassidis 5/3/95

if (nargin<6), flag=0; end

% Initialize Some Parameters and Pre-Allocate Space
[n,n1]=size(im);
qq=zeros(n,4);
i500=0;

% Main loop
for i=1:n,

% Display When Every 500th Point is Reached
if (i500==500), 
 disp(sprintf('      Program has reached point %5i',i-1))
 i500=0;
end
i500=i500+1;

% Get W and V Matrices
clear w v
[i1,j1]=find(av(i,:)==1);
for m=1:length(j1),
 w(:,m)=bm(i,j1(m)*3-2:j1(m)*3)'/norm(bm(i,j1(m)*3-2:j1(m)*3));
 v(:,m)=im(i,j1(m)*3-2:j1(m)*3)'/norm(im(i,j1(m)*3-2:j1(m)*3));
end

% Form Matrices for Eigenvector Decomposition
b=w*diag(wi(j1))*v';
s=b'+b;
z=[b(2,3)-b(3,2) b(3,1)-b(1,3) b(1,2)-b(2,1)];
sigg=trace(b);
k=[s-eye(3)*sigg z';z sigg];
[ve,e]=eig(k);
[hh,hh1]=max(diag(e));
qq(i,:)=ve(:,hh1)';

end

% Normalize Output
q1=qq(:,1);q2=qq(:,2);q3=qq(:,3);q4=qq(:,4); 
qnorm=(q1.*q1+q2.*q2+q3.*q3+q4.*q4).^0.5;
q=[q1./qnorm q2./qnorm q3./qnorm q4./qnorm];

% Try To Make Quaternions Continuous
if (flag==1),
 dq=abs(diff(q));[j,i]=max(max(dq));[i,j]=find(dq(:,i)>j-1.5);clear j
 si=size(i,1);l=length(q);if(round(si/2.01)*2~=si),i=[i;l];end;jj=length(i);
 for j=1:2:jj,q(i(j)+1:i(j+1),:)=-q(i(j)+1:i(j+1),:);end
end
