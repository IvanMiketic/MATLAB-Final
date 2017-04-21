%%%%%%%%%%%%%%%%%%%%%%%
%%%     MatLab      %%%
%%%   Final Exam    %%%
%%%    Tara Pena    %%%
%%%   Ivan Miketic  %%%
%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% These are some constants that we will be working with...

A = -1; % Original A which was given in the given potential
m = 1.419;  % Mass of an individual Rb atom in kg, techincally is multiplied by 10^-25
tspan = [1:25]; % time interval we will be using
run = 1000; run_A = [1:run]; % run is number of atoms, change one "run" if you wish to change # of atoms
a = 1; % This a is the radius of the deepest part of the well
radius = 1; radius2 = 1.5; radius3 = 2; % These three radii are for seeing the trajectories of the atoms with diff traps
k = 8.6173324e-5; T = 1; % Boltzmann constant, and T for temperature


%% Calculate trajectory without angular momentum

% This array will record the different numerical trajectories of the Rb
% atoms with different inital positions
rr = zeros(run, length(tspan));

% SwitchShella is a function which randomizes inital conditions and finds
% the trajectories of each atom and plots first 6 atom trajectories
[ pv, r, TT ] = SwitchShella(  tspan, A, radius, a, m, k, T,run, rr, 1);

% Using Maxwell Boltzmann probability, we can find the average position
% trajectory. The probabilites  are returned from the SwitchShellb function
figure(3)
plot(TT,pv)
caption = sprintf('Average Position with Radius %d', radius );
title(caption)
xlabel('t')
ylabel('<r>')

% We repeat the process with "radius2" to see what happens when the radius
% of well is changed.
rr2 = zeros(run, length(tspan));
[ pv1, r2, TT ] = SwitchShella(  tspan, A, radius2, a, m, k, T,run, rr2, 4);

figure(6)
plot(TT,pv1)
caption = sprintf('Average Position with Radius %d', radius2 );
title(caption)
xlabel('t')
ylabel('<r>')


% Again, the process is repeated, however with "radius3"
rr3 = zeros(run, length(tspan));
[ pv2, r3, TT ] = SwitchShella(  tspan, A, radius3, a, m, k, T,run, rr3, 7);

figure(9)
plot(TT,pv2)
caption = sprintf('Average Position with Radius %d', radius3 );
title(caption)
xlabel('t')
ylabel('<r>')

% We animated the trajectories of each atom with corresponding radii over
% time
figure(10)
filename = 'Rb_atoms_sphere_trap_nol.gif'; 
for i=1:length(TT)
    subplot(1,3,1)
   plot(TT.',r.')
   xlabel('time')
   ylabel('position')
   caption = sprintf('Trajectories - Radius of %d', radius );
   title(caption)
   text(0.1,190,sprintf('t = %1.2f',TT(i))); 
       subplot(1,3,2)
   plot(TT.',r2.')
   xlabel('time')
   ylabel('position')
   caption = sprintf('Radius of %d', radius2 );
   title(caption)
   text(0.1,190,sprintf('t = %1.2f',TT(i))); 
       subplot(1,3,3)
   plot(TT.',r3.')
   xlabel('time')
   ylabel('position')
   caption = sprintf('Radius of %d', radius3 );
   title(caption)
   text(0.1,190,sprintf('t = %1.2f',TT(i))); 
   frame = getframe(1);
   im = frame2im(frame); 
   [imind,cm] = rgb2ind(im,256);
   if i == 1;
       imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
   else
       imwrite(imind,cm,filename,'gif','WriteMode','append');
   end
end


%% Calculate trajectory with angular momentum 

% This array will record the different numerical trajectories of the Rb
% atoms with different inital positions
r = zeros(run,length(tspan));

% SwitchShella is a function which randomizes inital conditions and finds
% the trajectories of each atom and plots first 6 atom trajectories; the
% trajectories are different from SwitchShellb due to the use of angular
% momentum
[ pva, r1, TT ] = SwitchShellb( tspan, A, radius, a, m, k, T,run, r, 11);

% Using Maxwell Boltzmann probability, we can find the average position
% trajectory. The probabilites  are returned from the SwitchShellb function
figure(13)
plot(TT,pva)
caption = sprintf('Average Position with Radius %d', radius );
title(caption)
xlabel('t')
ylabel('<r>')

% We repeat the process with "radius2" to see what happens when the radius
% of well is changed.
r2 = zeros(run,length(tspan));
[ pvb, r2, TT ] = SwitchShellb( tspan, A, radius2, a, m, k, T,run, r2,14);

figure(16)
plot(TT,pvb)
caption = sprintf('Average Position with Radius %d', radius2 );
title(caption)
xlabel('t')
ylabel('<r>')

% Again, the process is repeated, however with "radius3"
r3 = zeros(run,length(tspan));
[ pvc, r3, TT ] = SwitchShellb( tspan, A, radius3, a, m, k, T,run, r3,17);

figure(19)
plot(TT,pvc)
caption = sprintf('Average Position with Radius %d', radius3 );
title(caption)
xlabel('t')
ylabel('<r>')

% We animated the trajectories of each atom with corresponding radii over
% time
figure(20)
filename = 'rb_atom_position_wl.gif'; 
for i=1:length(TT)
    subplot(1,3,1)
   plot(TT.',r1.')
   xlabel('time')
   ylabel('position')
   caption = sprintf('Trajectories Angular Momentum, Radius of %d', radius );
   title(caption)
   text(0.1,190,sprintf('t = %1.2f',TT(i))); 
       subplot(1,3,2)
   plot(TT.',r2.')
   xlabel('time')
   ylabel('position')
   caption = sprintf('Radius of %d', radius2 );
   title(caption)
      text(0.1,190,sprintf('t = %1.2f',TT(i)));
       subplot(1,3,3)
   plot(TT.',r3.')
   xlabel('time')
   ylabel('position')
   caption = sprintf('Radius of %d', radius3 );
   title(caption)
   text(0.1,190,sprintf('t = %1.2f',TT(i)));
   frame = getframe(1);
   im = frame2im(frame); 
   [imind,cm] = rgb2ind(im,256);
   if i == 1;
       imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
   else
       imwrite(imind,cm,filename,'gif','WriteMode','append');
   end
end

