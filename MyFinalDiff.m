function [ drdt ] = MyFinalDiff( t,r,A,m, radius, a )e

% MyFinalDiff sets up our differential equation 
% F = - Gradient ( V )
drdt = zeros(size(r));
x = r(1);
y = r(2);
z = r(3);
vx = r(4);
vy = r(5);
vz = r(6);

ra = sqrt(x.^2 + y .^2 + z.^2);

% This here is the  gradient of V in x,y,z
Fx = 4*A * x * ( ra - radius) ./ a.^2 * exp(-((ra - radius).^2)./(a.^2))./ ra.^2;
Fy = 4*A * y * ( ra - radius) ./ a.^2 * exp(-((ra - radius).^2)./(a.^2))./ ra.^2;
Fz = 4*A * z * ( ra - radius)./ a.^2 * exp(-((ra - radius).^2)./(a.^2))./ ra.^2;

drdt(1) = vx;
drdt(2) = vy;
drdt(3) = vz;

drdt(4) = Fx / m;
drdt(5) = Fy / m;
drdt(6) = Fz / m;
end
