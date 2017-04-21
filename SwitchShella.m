function [ splinepv, r , splinerangeb] = SwitchShella (tspan, A, radius, a, m, k,T, run, r, c)
% The for loop is for each atom trajectory
for i = 1:run
 % Randomize inital conditions for atom
[r_rand, theta_rand, phi_Rand, xo, yo, zo, vox, voy, voz, vo, ro] = randomizee(a, radius, k, T, m);
R0 = [xo, yo, zo, vox, voy, voz];
    
% ODE45 to solve for the trajectory
[TT,w] = ode45(@MyFinalDiff, tspan,R0,[],A,m, radius,a);

% Translate back into spherical .. w(:,1) = x, w(:,2) = y, w(:,3) = z
% Spline is to add accuracy to numerical values
RR = sqrt(w(:,1).^2 + w(:,2).^2 + w(:,3).^2);
splinerange = linspace(0,TT(end),length(TT));
splinedd = spline(TT,RR,splinerange);
theta =  atan(w(:,2)./RR)*180/pi;
splineee = spline(TT, theta, splinerange);
phi = acos(w(:,3)./RR);

% Calculate potential and total energy to use for Maxwell Boltzmann
U = -A * exp (-(r_rand - radius).^2/a^2);
E = .5 * m * vo.^2 + U;

% Probability of the trajectory of the atom
P_V = exp(-E.^2./2.*k.*T);
P_Vr(i) = exp(-E.^2./2.*k.*T) .* splinedd(1);

% Records atom's trajectory
r(i,:) = splinedd;

% Plots first six atom's trajectories in radial and angular directions
    if i <=6           
        figure(c)
        subplot(4,2,i)                       
        plot (splinerange, splinedd);                   
        caption = sprintf('Trial #%d of Radius %d',i, radius );   
        title(caption)                      
        xlabel('Time')
        ylabel('Position') 
        figure(c+1)
       subplot(6,1,i)                       
        plot (splinerange, splineee);                   
        caption = sprintf('Trial #%d of Radius %d',i, radius );  
        title(caption)                      
        xlabel('Time')
        ylabel('\Theta')
    end
end

% Calculate average trajectory
pv = P_Vr./sum(P_V);
splinerangeb = linspace(0, TT(end), length(pv));
splinepv = spline([1:run], pv, splinerangeb);
end
