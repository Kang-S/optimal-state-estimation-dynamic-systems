function [phi,theta,psi]=q2e(q,flag)
%function [phi,theta,psi]=q2e(q,flag)
%
% This m-file transforms quaternions to 
% Euler angles with the following rotations.
% Theta, phi, and psi are in radians.
%
%  The inputs are:
%    flag = 1 for 1-2-1
%         = 2 for 2-3-2
%         = 3 for 3-1-3
%         = 4 for 1-3-1
%         = 5 for 2-1-2
%         = 6 for 3-2-3
%         = 7 for 1-2-3
%         = 8 for 2-3-1
%         = 9 for 3-1-2
%         = 10 for 1-3-2
%         = 11 for 2-1-3
%         = 12 for 3-2-1

% Written by John L. Crassidis 4/24/95

q1=q(:,1);q2=q(:,2);q3=q(:,3);q4=q(:,4);
a=[q1.^2-q2.^2-q3.^2+q4.^2 2*(q1.*q2+q3.*q4) 2*(q1.*q3-q2.*q4) ...
2*(q1.*q2-q3.*q4) -q1.^2+q2.^2-q3.^2+q4.^2 2*(q2.*q3+q1.*q4) ... 
2*(q1.*q3+q2.*q4) 2*(q2.*q3-q1.*q4) -q1.^2-q2.^2+q3.^2+q4.^2];

if flag==1,
 t1=atan2(a(:,2),-a(:,3)) ;
 t2=acos(a(:,1));
 t3=atan2(a(:,4),a(:,7));
elseif flag==2,
 t1=atan2(a(:,6),-a(:,4));
 t2=acos(a(:,5));
 t3=atan2(a(:,8),a(:,2));
elseif flag==3,
 t1=atan2(a(:,7),-a(:,8));
 t2=acos(a(:,9));
 t3=atan2(a(:,3),a(:,6));
elseif flag==4,
 t1=atan2(a(:,3),a(:,2));
 t2=acos(a(:,1));
 t3=atan2(a(:,7),-a(:,4));
elseif flag==5,
 t1=atan2(a(:,4),a(:,6));
 t2=acos(a(:,5));
 t3=atan2(a(:,2),-a(:,8));
elseif flag==6,
 t1=atan2(a(:,8),a(:,7));
 t2=acos(a(:,9));
 t3=atan2(a(:,6),-a(:,3));
elseif flag==7,
 t1=atan2(-a(:,8),a(:,9));
 t2=asin(a(:,7));
 t3=atan2(-a(:,4),a(:,1));
elseif flag==8,
 t1=atan2(-a(:,3),a(:,1));
 t2=asin(a(:,2));
 t3=atan2(-a(:,8),a(:,5));
elseif flag==9,
 t1=atan2(-a(:,4),a(:,5));
 t2=asin(a(:,6));
 t3=atan2(-a(:,3),a(:,9));
elseif flag==10,
 t1=atan2(a(:,6),a(:,5));
 t2=asin(-a(:,4));
 t3=atan2(a(:,7),a(:,1));
elseif flag==11,
 t1=atan2(a(:,7),a(:,9));
 t2=asin(-a(:,8));
 t3=atan2(a(:,2),a(:,5));
elseif flag==12,
 t1=atan2(a(:,2),a(:,1));
 t2=asin(-a(:,3));
 t3=atan2(a(:,6),a(:,9));
end

phi=t1;theta=t2;psi=t3;
