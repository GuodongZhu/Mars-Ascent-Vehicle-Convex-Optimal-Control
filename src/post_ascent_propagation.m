%% Propagates equatorial orbit after optimal ascent phase has completed 
function xdot = post_ascent_propagation(x, G, M)
xdot = [x(4); %xdot1: x-velocity
        x(5); %xdot2: y-velocity
        x(6); %xdot3: z-velocity
        -(G*M*x(1))/((x(1)^2+x(2)^2+x(3)^2)^(3/2)); %xdot4: x-acceleration
        -(G*M*x(2))/((x(1)^2+x(2)^2+x(3)^2)^(3/2)); %xdot5: y-acceleration
        -(G*M*x(3))/((x(1)^2+x(2)^2+x(3)^2)^(3/2)); %xdot6: z-acceleration
        ];
end   
