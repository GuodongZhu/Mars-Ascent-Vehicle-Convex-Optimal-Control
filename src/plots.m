%% Producing plots of simulation output
function plots(t,r,v,m,u,R,D,r_orbit,G,M)

%Radial position vs. time
figure
plot(t,norms(r)-R)
grid on
xlabel('Time (s)')
ylabel('Radial Position (m)')
%print('-dpng','-r800','Images\radial_position')

%Total velocity vs. time
figure
plot(t,norms(v))
grid on
xlabel('Time (s)')
ylabel('Velocity Magnitude (m/s)')
%print('-dpng','-r800','Images\total_velocity')

%Radial velocity vs. time
figure
plot(t,dot(r,v)./norms(r))
grid on
xlabel('Time (s)')
ylabel('Radial Velocity (m/s)')
%print('-dpng','-r400','Images\radial_velocity')

%Mass vs. time
figure
plot(t,m)
grid on
xlabel('Time (s)')
ylabel('Mass (kg)')
%print('-dpng','-r400','Images\mass')

%Thrust components vs. time
figure
plot(t,u(1:3,:))
grid on
xlabel('Time (s)')
ylabel('Thrust Components (N)')
legend('T_x','T_y','T_z')
%print('-dpng','-r400','Images\thrust_components')

%Total thrust vs. time
figure
plot(t,u(4,:))
grid on
xlabel('Time (s)')
ylabel('Thrust Magnitude (N)')
%ylim([0 7500])
%print('-dpng','-r800','Images\total_thrust')

%Drag vs. time
figure
plot(t,D.*norms(v))
grid on
xlabel('Time (s)')
ylabel('Drag (N)')
%print('-dpng','-r400','Images\drag')

%% 3D Mars orbit
%Mars sphere and texture
figure
[X,Y,Z] = ellipsoid(0,0,0,R,R,R,40);
surf(X,Y,-Z)
mars_texture = imread('mars.jpg');
h = findobj('type','surface');
set(h,'cdata',mars_texture,'facecolor','texturemap','edgecolor','none')
hold on

%MAV optimal ascent trajectory 
plot3(r(1,:),r(2,:),r(3,:),'linewidth',2)
grid on
xlabel('x')
ylabel('y')
zlabel('z')

%Post-ascent equatorial orbit propagation 
x_init = [r(1,end), r(2,end), r(3,end), v(1,end), v(2,end), v(3,end)]
orbital_period = 2*pi*sqrt(r_orbit^3/(G*M));
t_span = [0 orbital_period];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,x] = ode15s(@(t,x) post_ascent_propagation(x, G, M), t_span, x_init, options);

plot3(x(:,1)',x(:,2)',x(:,3)','linewidth',2,'color',[0, 0.4470, 0.7410])
axis equal
hold off