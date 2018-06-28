function f=polfun_b(x,p,xf,c,k,q)

% Written by John L. Crassidis 9/03

% Function Routine for Van der Pol's Equation in Backwards Integration

fpart=[0 1;-4*c*xf(1)*xf(2)-k -2*c*(xf(1)^2-1)];
f=-[xf(2);-2*c*(xf(1)^2-1)*xf(2)-k*xf(1)]-[q*inv(p)+fpart]*([x(1);x(2)]-[xf(1);xf(2)]);