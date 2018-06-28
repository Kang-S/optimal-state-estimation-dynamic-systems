function mag_b=i2b(q,mag_i)
%function mag_b=i2b(q,mag_i)
%
% This m-file maps the inertial measurements to body
% measurements using the quaternion attitude matrix.
%
%  The input are:
%       q = quaternion (mx4)
%   mag_i = inertial measurement set (mx3)

% Written by John L. Crassidis 4/24/95

a1=q(:,1).^2-q(:,2).^2-q(:,3).^2+q(:,4).^2;
a2=2*(q(:,1).*q(:,2)+q(:,3).*q(:,4));
a3=2*(q(:,1).*q(:,3)-q(:,2).*q(:,4));
a4=2*(q(:,1).*q(:,2)-q(:,3).*q(:,4));
a5=-q(:,1).^2+q(:,2).^2-q(:,3).^2+q(:,4).^2;
a6=2*(q(:,2).*q(:,3)+q(:,1).*q(:,4));
a7=2*(q(:,1).*q(:,3)+q(:,2).*q(:,4));
a8=2*(q(:,2).*q(:,3)-q(:,1).*q(:,4));
a9=-q(:,1).^2-q(:,2).^2+q(:,3).^2+q(:,4).^2;

y1e=a1.*mag_i(:,1)+a2.*mag_i(:,2)+a3.*mag_i(:,3);
y2e=a4.*mag_i(:,1)+a5.*mag_i(:,2)+a6.*mag_i(:,3);
y3e=a7.*mag_i(:,1)+a8.*mag_i(:,2)+a9.*mag_i(:,3);

mag_b=[y1e y2e y3e];