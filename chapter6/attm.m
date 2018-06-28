function att=attm(q)
%function att=attm(q)
%
% This m-file determines the attitude matrix from a set of quaternions.
%
%  The input is:
%       q = quaternion (4x1)
%
%  The output as:
%     att = attitude matrix (3x3)

% John L. Crassidis 4/24/95

q=q(:)';
a1=q(:,1).^2-q(:,2).^2-q(:,3).^2+q(:,4).^2;
a2=2*(q(:,1).*q(:,2)+q(:,3).*q(:,4));
a3=2*(q(:,1).*q(:,3)-q(:,2).*q(:,4));
a4=2*(q(:,1).*q(:,2)-q(:,3).*q(:,4));
a5=-q(:,1).^2+q(:,2).^2-q(:,3).^2+q(:,4).^2;
a6=2*(q(:,2).*q(:,3)+q(:,1).*q(:,4));
a7=2*(q(:,1).*q(:,3)+q(:,2).*q(:,4));
a8=2*(q(:,2).*q(:,3)-q(:,1).*q(:,4));
a9=-q(:,1).^2-q(:,2).^2+q(:,3).^2+q(:,4).^2;
att=[a1 a2 a3;a4 a5 a6;a7 a8 a9];
