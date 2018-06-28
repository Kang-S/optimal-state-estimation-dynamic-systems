% This example shows the power of using a Kronecker product 
% least squares approach for a simple system. It outputs 
% the numerical accuracies of various approaches for 
% comparison purposes.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 1.10

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Gridded Points
for y = -2:0.2:2;
    rowy=[1 y y^2 y^3 y^4 y^5];
    if y == -2,
        hy=rowy;
    else
        hy=[hy;rowy];
    end
end
for x = -2:0.2:2;
    rowx=[1 x x^2 x^3 x^4 x^5];
    if x == -2,
        hx=rowx;
    else
        hx=[hx;rowx];
    end
end

% H Matrix and True Values
h=kron(hx,hy);
xtrue=ones(36,1);
ztrue=h*xtrue;

% Standard Least Squares
xhat1=inv(h'*h)*h'*ztrue;
norm_ls=norm(xhat1-xtrue)
disp(' ')

% % Pseudoinverse
% xhat2=pinv(h)*ztrue;
% norm_pseudo=norm(xhat2-xtrue)
% disp(' ')

% Kronecker Solution
xhat3=kron(inv(hx'*hx)*hx',inv(hy'*hy)*hy')*ztrue;
norm_kron=norm(xhat3-xtrue)
disp(' ')

% % Kronecker Solution with Pseudoinverse
% xhat4=kron(pinv(hx),pinv(hy))*ztrue;
% norm_kron_pseudo=norm(xhat4-xtrue)

% QR Solution
[q,r]=qr(h,0);
xhat5=inv(r)*q'*ztrue;
norm_qr=norm(xhat5-xtrue)
disp(' ')

% SVD Solution
[u,s,v]=svd(h,0);
xhat6=v*inv(s)*u'*ztrue;
norm_svd=norm(xhat6-xtrue)
disp(' ')
