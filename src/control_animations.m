%% Producing animations of simulation output
function control_animations(t_interp,r_ref,v_ref,m_ref,u_ff,r,v,m,u,R)

%Radial position vs. time
subplot(2,3,1)
plot(t_interp,norms(r_ref)-R)
grid on
xlabel('Time (s)')
ylabel('Radial Position (m)')
radial_position = animatedline('color', [0.8500, 0.3250, 0.0980]);

%Total velocity vs. time
subplot(2,3,2)
plot(t_interp,norms(v_ref))
grid on
xlabel('Time (s)')
ylabel('Total Velocity (m/s)')
total_velocity = animatedline('color', [0.8500, 0.3250, 0.0980]);

%Total thrust vs. time
subplot(2,3,6)
plot(t_interp,u_ff(4,:))
grid on
xlabel('Time (s)')
ylabel('Thrust Magnitude (N)')
total_thrust = animatedline('color', [0.8500, 0.3250, 0.0980]);

%Radial velocity vs. time
subplot(2,3,3)
plot(t_interp,dot(r_ref,v_ref)./norms(r_ref))
grid on
xlabel('Time (s)')
ylabel('Radial Velocity (m/s)')
radial_velocity = animatedline('color', [0.8500, 0.3250, 0.0980]);

%Mass vs. time
subplot(2,3,4)
plot(t_interp,m_ref)
grid on
xlabel('Time (s)')
ylabel('Mass (kg)')
mass = animatedline('color', [0.8500, 0.3250, 0.0980]);

%Control vs. time
subplot(2,3,5)
plot(t_interp,u_ff(1:3,:),'color',[0 0.4470 0.7410])
grid on
xlabel('Time (s)')
ylabel('Thrust Components (N)')
control_x = animatedline('color', [0.8500, 0.3250, 0.0980]);
control_y = animatedline('color', [0.8500, 0.3250, 0.0980]);
control_z = animatedline('color', [0.8500, 0.3250, 0.0980]);

%Animated figures
for k = 1:length(t_interp)
    addpoints(radial_position,t_interp(k),norms(r(:,k))-R)
    addpoints(total_velocity,t_interp(k),norms(v(:,k)))
    addpoints(radial_velocity,t_interp(k),dot(r(:,k),v(:,k))./norms(r(:,k)))
    addpoints(mass,t_interp(k),m(k))
    addpoints(control_x,t_interp(k),u(1,k))
    addpoints(control_y,t_interp(k),u(2,k))
    addpoints(control_z,t_interp(k),u(3,k))
    addpoints(total_thrust,t_interp(k),u(4,k))
    drawnow %limitrate
end