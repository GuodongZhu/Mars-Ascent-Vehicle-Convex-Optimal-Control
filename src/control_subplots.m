%% Producing subplots of simulation output
function control_subplots(t_interp,r_ref,v_ref,m_ref,u_ff,r,v,m,u,R)

%Radial position vs. time
subplot(2,3,1)
plot(t_interp,norms(r_ref)-R)
hold on
plot(t_interp,norms(r)-R)
hold off
grid on
xlabel('Time (s)')
ylabel('Radial Position (m)')
legend('r_{ref}','r','location','SE')
%print('-dpng','-r800','Images\radial_control')

%Total velocity vs. time
subplot(2,3,2)
plot(t_interp,norms(v_ref))
hold on
plot(t_interp,norms(v))
hold off
grid on
xlabel('Time (s)')
ylabel('Velocity Magnitude (m/s)')
legend('v_{ref}','v','location','SE')
%print('-dpng','-r800','Images\velocity_control')

%Total thrust vs. time
subplot(2,3,6)
plot(t_interp,u_ff(4,:))
hold on
plot(t_interp,u(4,:))
hold off
grid on
xlabel('Time (s)')
ylabel('Thrust Magnitude (N)')
legend('u_{ff}','u_{ff}+u_{fb}')
%ylim([0,7500])
%print('-dpng','-r800','Images\thrust_control')

%Radial velocity vs. time
subplot(2,3,3)
plot(t_interp,dot(r_ref,v_ref)./norms(r_ref))
hold on
plot(t_interp,dot(r,v)./norms(r))
hold off
grid on
xlabel('Time (s)')
ylabel('Radial Velocity (m/s)')

%Mass vs. time
subplot(2,3,4)
plot(t_interp,m_ref)
hold on
plot(t_interp,m)
hold off
grid on
xlabel('Time (s)')
ylabel('Mass (kg)')

%Control vs. time
subplot(2,3,5)
plot(t_interp,u_ff(1:3,:),'color',[0 0.4470 0.7410])
hold on
plot(t_interp,u(1:3,:),'color', [0.8500, 0.3250, 0.0980])
hold off
grid on
xlabel('Time (s)')
ylabel('Thrust Components (N)')