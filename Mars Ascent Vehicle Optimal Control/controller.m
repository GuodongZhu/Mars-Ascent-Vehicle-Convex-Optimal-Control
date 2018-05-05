%% Feedforward + PID feedback controller
function u = controller(r,r_ref,v,v_ref,m,u_ff,n,dt)
% previous_error = 0
% integral = 0
% loop:
%   error = setpoint - measured_value
%   integral = integral + error*dt
%   derivative = (error - previous_error)/dt
%   output = Kp*error + Ki*integral + Kd*derivative
%   previous_error = error
%   wait(dt)
%   goto loop

    %Using global variables to prevent them from being reset at each iteration
    global up_int_prev uv_int_prev e_prev %up_int uv_int
    
    %Gain scheduling is implemented by multiplying each gain by the MAV's current mass
    Kp = m*10; %Proportional gain
    %Ki = 0; %Integral gain
    Kd = 0; %Integral gain
    Ki = m*20; %Integral gain
    %Kd = m*0.5; %Derivative gain
%     Kp = m*0.01; %Proportional gain
%     Ki = m*0.0; %Integral gain
%     Kd = m*0.01; %Derivative gain
    
    %Position and velocity errors
    e = [r_ref-r; v_ref-v];
    
    if n == 1
        up_int_prev = zeros(3,1);
        uv_int_prev = zeros(3,1);
        e_prev = zeros(6,1);
    end
    
    %Proportional control
    up_prop = Kp*e(1:3);
    uv_prop = Kp*e(4:6);
    %Integral control
    up_int = up_int_prev + Ki*e(1:3)*dt;
    uv_int = uv_int_prev + Ki*e(4:6)*dt;
    %Derivative control
    up_der = Kd*(e(1:3) - e_prev(1:3))/dt;
    uv_der = Kd*(e(4:6) - e_prev(4:6))/dt;
    
    %PID feedback control
    u_fb = up_prop + uv_prop + up_der + uv_der + up_int + uv_int;
      
    %Combined control using optimal feedforward and feedback control inputs
    u(1:3,1) = u_ff(1:3,1) + u_fb;
    u(4,1) = norms(u(1:3,1));
    
    up_int_prev = up_int;
    uv_int_prev = uv_int;
    e_prev = e;
end