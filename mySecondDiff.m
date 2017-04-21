function drdt = mySecondDiff( t, r, A,m, radius,a,L)

% mySecondDiff sets up our differential equation
% F = - Gradient ( V )
drdt = zeros(size(r));
x = r(1);
y = r(2);
z = r(3);
vx = r(4);
vy = r(5);
vz = r(6);

ra = sqrt(x.^2 + y.^2 + z.^2);

% Very similar to MyFinalDiff, however we have the L^2/mr^2 term for
% angular momentum.
Fx = 4*A * x * ( ra - radius) ./ a.^2 * exp(-((ra - radius).^2)./(a.^2))./ ra.^2 - x *L.^2./2*m.^2./ra.^4;
Fy = 4*A * x * ( ra - radius) ./ a.^2 * exp(-((ra - radius).^2)./(a.^2))./ ra.^2 - y *L.^2./2*m.^2./ra.^4;
Fz = 4*A * x * ( ra - radius) ./ a.^2 * exp(-((ra - radius).^2)./(a.^2))./ ra.^2 - z *L.^2./2*m.^2./ra.^4;

drdt(1) = vx;
drdt(2) = vy;
drdt(3) = vz;

drdt(4) = Fx / m;
drdt(5) = Fy / m;
drdt(6) = Fz / m;
end