hold off
clf
clc
clear

for n = [1,2,3]
    a = 0.1; % Boundary one
    b = 10; % infinity [1,2,3]
    ph = 4; %bisection upper bound
    pl = -4; %Bisection lower bound
    %p = (ph + pl)/2; %bisection first mid point
    p = 0.2; %Arbitray initial shooting parameter used for NM
    tol = 1; %Set initial tolerance
    i = 1; %Counter for number of iterations
    u0 = [0,p,0,1]; %Set initial conditions for our system
    step_size = 0.2; %Set step size 
    while i < 1000 && tol > 1e-7
        %Integrate the system using runge-kutta method
        [r,u] = ode45(@(r,u) shoot(r,u,n), [a:step_size:b], u0,n);
        tol = abs(1-u(end,1)) %Display the tolerance
        tor = u(end,1); %Purely aesthetic
        %Bisection scheme for updating p.
        if(tor > 1)
            ph = p; %Move upper bound
            %p = (pl + ph)/2; %Bisect
            %u0 = [0,p,0,1]; %Update IC's with new p value
        end
        if(tor < 1) %Note that we exclude tor = 1 explicitly as if tor = 1 we are finished.
            pl = p; %Move lower bound
            %p = (pl + ph)/2; %Bisect
            %u0 = [0,p,0,1]; %Update IC's with new p value
        end
        %End of Bisection scheme
        
        %NM scheme. Note that we use bisection primarily for stability
        p = p - ((u(end,1)-1)/(u(end,3))); %Update p according to NM scheme
        u0 = [0,p,0,1]; %Update the IC's with new p value
        
        i = i + 1; %Update iteration counter
        figure(1);plot(r,u(:,1)); %Plot the updated solution curve
        hold off;
        pause(0.1);
    end
    fig2 = figure(2);
    plot(r,u(:,1))
    hold on;
    legend("n=1", "n=2", "n=3");
    title("Nonlinear shooting for a = 0.1, b = 15, tol = 1e-7");
    xlabel('r');
    ylabel('u');
    pause(1);
end

%The system of equations for u and z. Note z is only needed for NM scheme
function [du] = shoot(r,u, n)
    uu = u(1);
    v = u(2);
    z = u(3);%Extract the IC's
    w = u(4);
    %Update scheme needed for integrating u
    du(1,1) = v;
    du(2,1) = (-1/r)*v  - (uu/(1-uu^2))*(v^2-(n^2)/(r^2)) - uu*(1-uu^2);
    %Update scheme needed for integrating z
    du(3,1) = w;
    du(4,1) = z*( -((1+uu^2)/((1-uu^2)^2))*(v^2 - (n^2)/(r^2)) - 1 + 3*uu^2) + w*(-1/r  - ((2*uu)*v)/(1-uu^2));
end


