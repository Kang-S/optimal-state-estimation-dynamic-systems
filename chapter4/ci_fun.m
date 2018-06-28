function f=ci_fun(omega,p1,p2)

f=trace(inv(omega*inv(p1)+(1-omega)*inv(p2)));