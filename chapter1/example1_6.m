% This example shows Newton's root solver to determine 
% the roots of a third-order polynomial. It presents the 
% iteration results for three different starting values.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 1.6

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Starting Value of 0
xit=zeros(8,1);xit(1)=0;
for i = 1:7
 f=xit(i)^3+6*xit(i)^2+11*xit(i)+6;
 dfdx=3*xit(i)^2+12*xit(i)+11;
 xit(i+1)=xit(i)-inv(dfdx)*f;
end
iteration_results=xit

disp(' Press any key to continue')
disp(' ')
pause

% Starting Value of -1.6
xit=zeros(8,1);xit(1)=-1.6;
for i = 1:7
 f=xit(i)^3+6*xit(i)^2+11*xit(i)+6;
 dfdx=3*xit(i)^2+12*xit(i)+11;
 xit(i+1)=xit(i)-inv(dfdx)*f;
end
iteration_results=xit

disp(' Press any key to continue')
disp(' ')
pause

% Starting Value of -5.0
xit=zeros(8,1);xit(1)=-5.0;
for i = 1:7
 f=xit(i)^3+6*xit(i)^2+11*xit(i)+6;
 dfdx=3*xit(i)^2+12*xit(i)+11;
 xit(i+1)=xit(i)-inv(dfdx)*f;
end
iteration_results=xit