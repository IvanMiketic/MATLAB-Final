function [r_rand, theta_rand, phi_Rand, xo, yo, zo, vox, voy, voz, vo, ro] = randomizee( a, radius, k, T, m)
% The shape of our trap is a sphere
% Therefore we have these rand functions for our random inital conditions
r_rand = (radius - a) + 2*a*rand;
theta_rand = pi*rand;
phi_Rand = 2*pi*rand;

% Conversions into cartesian
xo = r_rand*sin(phi_Rand)*cos(theta_rand);
yo = r_rand*sin(phi_Rand)*sin(theta_rand);
zo = r_rand*cos(phi_Rand);

% Find the max velocity due to mass and temperature
v_max = 2*sqrt(k*T/m);
vox = (rand-0.5)*2*v_max;voy = (rand-0.5)*2*v_max; voz = (rand-0.5)*2*v_max;
vo = sqrt(vox.^2 + voy.^2 + voz.^2);
ro = sqrt(xo.^2 + yo.^2 + zo.^2);
end