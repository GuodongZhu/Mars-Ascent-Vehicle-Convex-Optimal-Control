% -------------------------------------------------------------------------
% The techniques used within this algorithm were based on similar methods
% implemented in the PhD thesis by Xinfu Liu titled "Autonomous Trajectory
% Planning by Convex Optimization", Iowa State University, 2013.
% 
% The following code implements a convex optimization approach to allow
% the Mars Ascent Vehicle (MAV) to autonomously plan and execute a fuel
% optimal ascent trajectory, from the surface of Mars to a stable circular
% orbit.
% 
% The optimal guidance algorithm makes use of convex relaxations and
% sequential convex programming (also known as successive convexification)
% to solve a series of convex subproblems. This enables rapid progess
% towards an optimal solution that minimizes fuel use along the MAV's
% ascent trajectory.
% 
% Further techniques were also implemented from the following papers:
% 
% 1) Solving Nonconvex Optimal Control Problems by Convex Optimization. Liu
% X., Lu P., 2013
% 
% 2) Successive Convexification for Fuel-Optimal Powered Landing with
% Aerodynamic Drag and Non-Convex Constraints. Szmuk M., Acikmese B.,
% Berning A., 2016
% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
% Note: Algorithm requires the following software:
% 1) YALMIP optimization environment (https://yalmip.github.io/)
% 2) MOSEK optimization package (https://www.mosek.com/)
% 3) A norms function to allow for computation of multiple vector norms
% -------------------------------------------------------------------------

clear
clc
close all
format long g

%% Simulation parameters
k_max = 20; %Maximum number of SOCP solution iterations
tf = 800; %Simulation end-time
N = 100; %Number of time nodes from t0 to tf
dt = tf/N; %Time step-size between each node

G = 6.672e-11; %Gravitational constant
M = 6.39e23; %Mass of Mars (kg)
R = 3.39e6; %Radius of Mars (m)

GLOM = 346; %Gross lift-off mass (kg)
m_prop = 259.7; %Propellant mass of stage 1 (kg)
m_payload = 18; %Payload mass (kg)
m_struct = GLOM - m_prop - m_payload; %Structural mass of stage 1 (kg) (includes payload of 5kg)
Isp = 314; %Engine specific impulse (s)
g0 = 9.81; %Reference acceleration (m/s^2)
%T_max = 7000; %Maximum thrust level (N)
T_max = 10000; %Maximum thrust level (N)

r_orbit = R+479000; %Desired orbital altitude above Mars' surface (m)
v_orbit = sqrt((G*M)/r_orbit); %Orbital velocity required to maintain circular orbit at altitude r = rorbit

%Initial position in spherical coordinates
%Coordinate convention based on http://mathworld.wolfram.com/SphericalCoordinates.html
phi0 = 80; %phi = 80 => latitude = 10 degrees
theta0 = 90;
r0 = R;

%Spherical to Cartesian coordinate transformation
x0 = r0*cosd(theta0)*sind(phi0) %Initial x-position (m)
y0 = r0*sind(theta0)*sind(phi0) %Initial y-position (m)
z0 = r0*cosd(phi0) %Initial z-position (m)

vx0 = 0; %Initial x-velocity (m/s)
vy0 = 0; %Initial y-velocity (m/s)
vz0 = 0; %Initial z-velocity (m/s)

%Initial conditions
r0 = [x0; y0; z0]; %Initial position (m)
v0 = [vx0; vy0; vz0]; %Initial velocity (m/s)
m0 = GLOM; %Initial mass (kg)
mc0 = log(m0); %Initial transformed mass

%Final conditions
xf = (r_orbit)/5; %Final x-position (m), obtained from Mars gravity turn burn-coast-burn ascent
yf = R+(r_orbit)/5; %Final y-position (m), obtained from Mars gravity turn burn-coast-burn ascent
zf = 0; %Final z-position (m)
vxf = v_orbit/2; %Final x-velocity (m/s)
vyf = -v_orbit/2; %Final y-velocity (m/s)
vzf = 0; %Final z-velocity (m/s)
rf = [xf; yf; zf]; %Final position (m)
vf = [vxf; vyf; vzf]; %Final velocity (m/s)
mf = m_struct + m_payload; %Finall mass (kg)

%Constant B matrix
v_ex = g0*Isp;
B = [
    0 0 0 0;
    0 0 0 0;
    0 0 0 0;
    1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 -1/v_ex];

%Initial radial position vector, r, at iteration k = 1 FOR EACH time node
r1 = [linspace(r0(1),rf(1),N); linspace(r0(2),rf(2),N); linspace(r0(3),rf(3),N)];

%Initial velocity vector, v, at iteration k = 1 FOR EACH time node
v1 = [linspace(v0(1),vf(1),N); linspace(v0(2),vf(2),N); linspace(v0(3),vf(3),N)];

%Initial mass vector approximation, z, at iteration k = 1 FOR EACH time node n
mc1 = linspace(log(m0), log(mf), N);

%Initial drag vector, D, at iteration k = 1 FOR EACH time node
S_D = pi*(0.57/2)^2; %Drag reference area
C_D = 0.295; %Coefficient of drag

%Initial A matrix at iteration k = 1 FOR EACH time node
Ak = zeros(7,7,N);
D1 = drag(r1,v1,R,S_D,C_D);
m1 = exp(mc1);
for t = 1:N
    Ak(:,:,t) = [
    0 0 0 1 0 0 0;
    0 0 0 0 1 0 0;
    0 0 0 0 0 1 0;
    -(G*M)/(norm(r1(:,t)).^3) 0 0 -D1(t)/m1(t) 0 0 0;
    0 -(G*M)/(norm(r1(:,t)).^3) 0 0 -D1(t)/m1(t) 0 0;
    0 0 -(G*M)/(norm(r1(:,t)).^3) 0 0 -D1(t)/m1(t) 0;
    0 0 0 0 0 0 0];
end

rk = zeros(1,N); %Arbitrary initial condition for r^(k-1)(t)
vk = zeros(1,N); %Arbitrary initial condition for v^(k-1)(t)
mck = zeros(1,N); %Arbitrary initial condition for mc^(k-1)(t)
yk = [r1; v1; mc1]; %Arbitrary initial condition for y^(k-1)(t)
epsilon = (1e-1)*ones(7,N); %Stopping criteria for state vector y
delta = [(1e3)*ones(1,N); (1e3)*ones(1,N); (1e3)*ones(1,N); (1e2)*ones(1,N); (1e2)*ones(1,N); (1e2)*ones(1,N); (1e3)*ones(1,N)]; %Trust-region radius
gamma = 0.01;

%% Scaling factors
length_scale = 1/r_orbit;
velocity_scale = 1/v_orbit;
mass_scale = 1/mc0;
time_scale = length_scale/velocity_scale;
acceleration_scale = m0/T_max;
        
rx_scale = 1/xf;
ry_scale = 1/yf;
rz_scale = 1;
vx_scale = 1/vxf;
vy_scale = 1/vyf;
vz_scale = 1;
vr_scale = 1/1370;

r_scale = [rx_scale 0 0
            0 ry_scale 0
            0 0 rz_scale];
        
v_scale = [vx_scale 0 0
              0 vy_scale 0
              0 0 vz_scale];
          
M_scale = [rx_scale 0 0 0 0 0 0;
            0 ry_scale 0 0 0 0 0;
            0 0 rz_scale 0 0 0 0;
            0 0 0 vx_scale 0 0 0;
            0 0 0 0 vy_scale 0 0;
            0 0 0 0 0 vz_scale 0;
            0 0 0 0 0 0 mass_scale];

objective_value = [];
trajectory_convergence_history = [];
mass_convergence_history = [];
velocity_convergence_history = [];
control_convergence_history = [];

%% Optimization algorithm
tic
options = sdpsettings('solver','mosek','verbose',1)
for k = 1:k_max
    k
	%Defining optimization variables
    y = sdpvar(7,N);
    u = sdpvar(4,N);
    
    constraints = [];

    constraints = [
        %Initial conditions
        M_scale*y(:,1) == M_scale*[r0; v0; mc0]
        ];
        
        %Linearized terminal constraints
        if k == 1
        constraints = [
            constraints;
            
            %Linearized position constraint
            length_scale*(norm(r1(:,N)) - r_orbit + (r1(1,N)*(y(1,N) - r1(1,N)))/norm(r1(:,N)) + (r1(2,N)*(y(2,N) - r1(2,N)))/norm(r1(:,N)) + (r1(3,N)*(y(3,N) - r1(3,N)))/norm(r1(:,N))) == 0;%gamma*((r1(1,N)^2 + r1(2,N)^2 + r1(3,N)^2)^(1/2) - r_orbit)
            
            %Linearized velocity constraint
            velocity_scale*(norm(v1(:,N)) - v_orbit + (v1(1,N)*(y(4,N) - v1(1,N)))/norm(v1(:,N)) + (v1(2,N)*(y(5,N) - v1(2,N)))/norm(v1(:,N)) + (v1(3,N)*(y(6,N) - v1(3,N)))/norm(v1(:,N))) == 0;%gamma*((v1(1,N)^2 + v1(2,N)^2 + v1(3,N)^2)^(1/2) - v_orbit)
            
            %Linearized flight-path angle constraint
            vr_scale*(r1(1,N)*v1(1,N) + r1(2,N)*v1(2,N) + r1(3,N)*v1(3,N) + v1(1,N)*(y(1,N) - r1(1,N)) + v1(2,N)*(y(2,N) - r1(2,N)) + v1(3,N)*(y(3,N) - r1(3,N)) + r1(1,N)*(y(4,N) - v1(1,N)) + r1(2,N)*(y(5,N) - v1(2,N)) + r1(3,N)*(y(6,N) - v1(3,N))) == 0;%gamma*(v1(1,N)*r1(1,N) + v1(2,N)*r1(2,N) + v1(3,N)*r1(3,N))
            
            %Linearized orbital inclination constraint
            %(r1(1,N)*v1(2,N) - r1(2,N)*v1(1,N)) + v1(2,N)*(r(1,N) - r1(1,N)) - v1(1,N)*(r(2,N) - r1(2,N)) - r1(2,N)*(v(1,N) - v1(1,N)) + r1(1,N)*(v(2,N) - v1(2,N)) - r_orbit*v_orbit == 0;%cosd(i_orbit) == 0;
            rz_scale*r1(3,N) == rz_scale*0;
            vz_scale*v1(3,N) == vz_scale*0;
            ];
        else
        constraints = [
            constraints;
            
            %Linearized position constraint
            length_scale*(norm(rk(:,N)) - r_orbit + (rk(1,N)*(y(1,N) - rk(1,N)))/norm(rk(:,N)) + (rk(2,N)*(y(2,N) - rk(2,N)))/norm(rk(:,N)) + (rk(3,N)*(y(3,N) - rk(3,N)))/norm(rk(:,N))) == length_scale*0;%gamma*((r1(1,N)^2 + r1(2,N)^2)^(1/2) - r_orbit)
            
            %Linearized velocity constraint
            velocity_scale*(norm(vk(:,N)) - v_orbit + (vk(1,N)*(y(4,N) - vk(1,N)))/norm(vk(:,N)) + (vk(2,N)*(y(5,N) - vk(2,N)))/norm(vk(:,N)) + (vk(3,N)*(y(6,N) - vk(3,N)))/norm(vk(:,N))) == velocity_scale*0;%gamma*((v1(1,N)^2 + v1(2,N)^2)^(1/2) - v_orbit)
            
            %Linearized flight-path angle constraint
            vr_scale*(rk(1,N)*vk(1,N) + rk(2,N)*vk(2,N) + rk(3,N)*vk(3,N) + vk(1,N)*(y(1,N) - rk(1,N)) + vk(2,N)*(y(2,N) - rk(2,N)) + vk(3,N)*(y(3,N) - rk(3,N)) + rk(1,N)*(y(4,N) - vk(1,N)) + rk(2,N)*(y(5,N) - vk(2,N)) + rk(3,N)*(y(6,N) - vk(3,N))) == velocity_scale*0;%gamma*(v1(1,N)*r1(1,N) + v1(2,N)*r1(2,N) + v1(3,N)*r1(3,N))
            
            %Linearized orbital inclination constraint
            %(rk(1,N)*vk(2,N) - rk(2,N)*vk(1,N)) + vk(2,N)*(r(1,N) - rk(1,N)) - vk(1,N)*(r(2,N) - rk(2,N)) - rk(2,N)*(v(1,N) - vk(1,N)) + rk(1,N)*(v(2,N) - vk(2,N)) - hf == 0;%*cosd(i_orbit) == 0;
            rz_scale*y(3,N) == rz_scale*0; 
            vz_scale*y(6,N) == vz_scale*0; 
            ];  
        end

    %Dynamic equality constraints at each time node N
    for t = 1:N-1
            
        %Runge-Kutta 4th-order coefficients
        k1 = Ak(:,:,t)*y(:,t) + B*u(:,t);
        k2 = Ak(:,:,t)*(y(:,t) + (dt*k1)/2) + B*u(:,t);
        k3 = Ak(:,:,t)*(y(:,t) + (dt*k2)/2) + B*u(:,t);
        k4 = Ak(:,:,t)*(y(:,t) + dt*k3) + B*u(:,t);

    constraints = [
        constraints;
        %Runge-Kutta 4th-order integration scheme
        M_scale*y(:,t+1) == M_scale*(y(:,t) + (dt/6)*(k1+ 2*k2 + 2*k3 + k4));
        ];
    end
    
    %Static equality and inequality constraints at each node N
    constraints = [
        constraints;
        acceleration_scale*u(4,:) >= acceleration_scale*0;
        acceleration_scale*u(4,:) <= acceleration_scale*(3*g0); %Maximum acceleration constraint
        mass_scale*y(7,:) <= mass_scale*mc0; %Ensures not all fuel is used, seems to produce error if constraint included when k = 1            
        cone([acceleration_scale*u(4,:); acceleration_scale*u(1:3,:)]); 
        ];

    if k == 1
        constraints = [
            constraints;
            acceleration_scale*u(4,:) <= acceleration_scale*(T_max.*exp(-mc1).*(1-(y(7,:)-mc1)));
            ];
    else
        constraints = [
            constraints;
            acceleration_scale*u(4,:) <= acceleration_scale*(T_max.*exp(-mck).*(1-(y(7,:)-mck)));
            ];
    end
    
    %Performing optimization
    objective = mass_scale*(-y(7,N)); %Maximizing final mass
    optimize(constraints,objective,options);

    y_opt = y;
    u_opt = u;
    
    r = value(y_opt(1:3,:));
    v = value(y_opt(4:6,:));
    mc = value(y_opt(7,:));
    y = value(y_opt);
    u = value(u_opt);
    
    %Position testing
    t = linspace(0,tf,N);
    hold on
    plot(t,norms(r1)-R,'k')
    plot(t,norms(r)-R)
    xlabel('Time (s)')
    ylabel('Radial Position (m)')
    grid on
    ylim([0,500000])
    hold off
    
    %Mars atmospheric drag model for iteration k+1
    D = drag(r,v,R,S_D,C_D);
    
    %Computing new constant A matrix from iteration k to be used in iteration k+1  
    m = exp(mc);
    for t = 1:N
        Ak(:,:,t) = [
        0 0 0 1 0 0 0;
        0 0 0 0 1 0 0;
        0 0 0 0 0 1 0;
        -(G*M)/(norm(r(:,t)).^3) 0 0 -D(t)/m(t) 0 0 0;
        0 -(G*M)/(norm(r(:,t)).^3) 0 0 -D(t)/m(t) 0 0;
        0 0 -(G*M)/(norm(r(:,t)).^3) 0 0 -D(t)/m(t) 0;
        0 0 0 0 0 0 0];
    end
   
    delta_y = abs(y-yk);
    if k == 1
        delta_m = abs(exp(mc(N))-exp(mc1(N)));
    else 
        delta_m = abs(exp(mc(N))-exp(mck(N)));
    end
    
    %Recording results for convergence plotting
    objective_value = [objective_value, value(objective)]
    trajectory_convergence_history = [trajectory_convergence_history; norms(r)];
    mass_convergence_history = [mass_convergence_history; m];
    velocity_convergence_history = [velocity_convergence_history; norms(v)];
    control_convergence_history = [control_convergence_history; norms(u(1:3,:))];
    
    %Stopping criteria: loop exits if change between solutions k and k-1 is small enough
    %Boundary constraints must also be met
    if delta_m < 1 && abs(norm(y(1:3,end)) - r_orbit) < 0.1 && abs(norm(y(4:6,end)) - v_orbit) < 0.1
        break
    else
        yk = y; %Computing y^(k-1)(t)
        rk = r; %Computing r^(k-1)(t)
        vk = v; %Computing v^(k-1)(t)
        mck = mc; %Computing mc^(k-1)(t)
    end
end
toc

%% Plotting convergence results
t = linspace(0,tf,N);

%Position testing
figure
hold on
plot(t,norms(r1)-R,'k')
plot(t,trajectory_convergence_history-R)
xlabel('Time (s)')
ylabel('Radial Position (m)')
grid on
ylim([0,500000])
hold off
legend('Initial guess','Iteration 1', 'Iteration 2', 'Iteration 3', 'Iteration 4', 'Iteration 5', 'Iteration 6', 'Iteration 7', 'Location', 'NW')
%print('-dpng','-r800','Images\trajectory_convergence_1')

%Velocity convergence history
figure
hold on
plot(t,norms(v1),'k')
plot(t,velocity_convergence_history)
xlabel('Time (s)')
ylabel('Total Velocity (m/s)')
grid on
legend('Initial guess','Iteration 1','Iteration 2', 'Iteration 3', 'Iteration 4', 'Iteration 5', 'Iteration 6', 'Iteration 7', 'Location', 'Best')
hold off
%print('-dpng','-r800','Images\total_velocity_convergence')

%Mass convergence history
figure
hold on
plot(t,exp(mc1),'k')
plot(t,mass_convergence_history)
xlabel('Time (s)')
ylabel('Mass (kg)')
grid on
legend('Initial guess','Iteration 1','Iteration 2', 'Iteration 3', 'Iteration 4', 'Iteration 5', 'Iteration 6', 'Iteration 7', 'Location', 'Best')
hold off
%print('-dpng','-r800','Images\mass_convergence')

%Control convergence history
figure
%hold on
plot(t,control_convergence_history.*exp(log(mass_convergence_history)))
xlabel('Time (s)')
ylabel('Thrust Magnitude (N)')
grid on
legend('Iteration 1','Iteration 2', 'Iteration 3', 'Iteration 4', 'Iteration 5', 'Iteration 6', 'Iteration 7', 'Location', 'Best')
%ylim([0, 7500])
%print('-dpng','-r800','Images\thrust_convergence')

%Objective value
figure
plot(objective_value,'o-')
xlabel('Solution Iteration')
ylabel('Objective Value')
grid on
%print('-dpng','-r800','Images\objective_convergence')

%% Plotting output
T = zeros(4,N);
for i = 1:N
    T(:,i) = u(:,i)*m(i);
end

plots(t,r,v,m,T,R,D, r_orbit,G,M)

%% Feedforward + feedback control implementation
%Computing splines of generated trajectory for use by controller
t = linspace(0,tf,N);
r_spline = spapi(4, t, r);
v_spline = spapi(3, t, v);
u_spline = spapi(2, t, u);
m_spline = spapi(2, t, m);

%Time vector used for spline plots
N_interp = 100*N;
t_interp = linspace(0,tf,N_interp);

%Reference spline plots
r_ref = fnval(r_spline, t_interp); %Desired reference position
v_ref = fnval(v_spline, t_interp); %Desired reference velocity
m_ref = fnval(m_spline, t_interp);
u_ref = fnval(u_spline, t_interp); %Feedforward control input

%Feedforward thrust control input
for i = 1:length(m_ref)
    T_ref(:,i) = u_ref(:,i)*m_ref(i);
end

%% Running simulation using optimized control inputs
%Clearing out old states and controls from optimization above
y = [];
r = [];
v = [];
u = [];
k1 = [];
k2 = [];
k3 = [];
k4 = [];

%Initial conditions
y = [x0; y0; z0; vx0; vy0; vz0; m0];
r = [x0; y0; z0];
v = [vx0; vy0; vz0];
dt = tf/N_interp; %Time step-size between each node
u_ff = T_ref; %Feedforward thrust control input
u = zeros(4,length(u_ref));
%C_D = 2 %Testing model uncertainty with different C_D value
for n = 1:length(t_interp)-1
    %Computing control input at each iteration n
    u(:,n) = controller(r(:,n),r_ref(:,n),v(:,n),v_ref(:,n),m(n),u_ff(:,n),n,dt);

    %Computing drag components
    D = drag(r,v,R,S_D,C_D);
    
    %Variable state A matrix
    %m(n) = exp(m(n));
    A(:,:,n) = [
    0 0 0 1 0 0 0;
    0 0 0 0 1 0 0;
    0 0 0 0 0 1 0;
    -(G*M)/(norm(r(:,n)).^3) 0 0 -D(n)/m(n) 0 0 0;
    0 -(G*M)/(norm(r(:,n)).^3) 0 0 -D(n)/m(n) 0 0;
    0 0 -(G*M)/(norm(r(:,n)).^3) 0 0 -D(n)/m(n) 0;
    0 0 0 0 0 0 0];
    
    %Constant control B matrix
    B = [
        0 0 0 0;
        0 0 0 0;
        0 0 0 0;
        1/m(n) 0 0 0;
        0 1/m(n) 0 0;
        0 0 1/m(n) 0;
        0 0 0 -1/v_ex];
    
    %Runge-Kutta 4th-order integration scheme
    k1 = A(:,:,n)*y(:,n) + B*u(:,n);
    k2 = A(:,:,n)*(y(:,n) + (dt*k1)/2) + B*u(:,n);
    k3 = A(:,:,n)*(y(:,n) + (dt*k2)/2) + B*u(:,n);
    k4 = A(:,:,n)*(y(:,n) + dt*k3) + B*u(:,n);
    y(:,n+1) = y(:,n) + (dt/6)*(k1+ 2*k2 + 2*k3 + k4);
                
    r(:,n+1) = y(1:3,n+1);
    v(:,n+1) = y(4:6,n+1);
    m(n+1) = y(7,n+1);
    
    %Filling in final column of drag and control matrices/vectors
    if n == length(t_interp)-1
        D = drag(r,v,R,S_D,C_D);
        u(:,n+1) = zeros(4,1);
    end
end

%% Error between actual and reference trajectories
figure
plot(t_interp,norms(r_ref)-norms(r))
grid on
xlabel('Time (s)')
ylabel('Position Error (m)')

figure
plot(t_interp,norms(v_ref)-norms(v))
grid on
xlabel('Time (s)')
ylabel('Velocity Error (m)')

figure
plot(t_interp,m_ref-m)
grid on
xlabel('Time (s)')
ylabel('Mass Error (kg)')

%% Subplots of output trajectory
control_subplots(t_interp,r_ref,v_ref,m_ref,u_ff,r,v,m,u,R)

%% Animations of output trajectory
%control_animations(t_interp,r_ref,v_ref,m_ref,u_ff,r,v,m,u,R)