hold off
clf
clc
clear
hold on



file = fopen("inf.txt", 'w');
for n = [1,2,3]
     fprintf(file,'%6.2f %12.8f\n', n);
    clearvars -except n file
for step_size = [0.002:0.001:0.5] %Variable to vary
    a = 0.1; % Boundary one
    b = 5; % infinity
    ph = 4; %bisection
    pl = -4;
    p = (ph + pl)/2; %bisection
    %p = 0.2; %Arbitr0.001ay initial shooting parameter
    tol = 1;
    i = 1;
    u0 = [0,p,0,1];
while i < 100 && tol > 1e-7
   [r,u] = ode45(@(r,u) shoot(r,u,n), [a:step_size:b], u0,n);
   tol = abs(1-u(end,1));
   tor = u(end,1);
   if(tor > 1)
      ph = p;
      p = (pl + ph)/2;
      u0 = [0,p,0,1];   
   end
   if(tor < 1)
      pl = p;
      p = (pl + ph)/2;
      u0 = [0,p,0,1];  
   end
   
   %p = p - ((u(end,1)-1)/(u(end,3)));
   %u0 = [0,p,0,1];

   i = i + 1;
end %end of while loop
    bool = ~isnan(tol) && tol < 1e-7;%tolerence to be considered stable
    fprintf(file,'%f %f\n', step_size, bool);
end  %end of vrbl loop 

   
    %plot(r,u(:,1))
    %legend("n=1", "n=2", "n=3")
    %title("Nonlinear shooting for a = 0.1, b = 15, tol = 1e-7")
    %xlabel('r')
    %ylabel('u')
end %end of n loop
fprintf(file, "TESTING B OVER");
fclose(file);
function [du] = shoot(r,u, n)
    %n = 2;
    uu = u(1);
    v = u(2);
    z = u(3);
    w = u(4);
    
    
    du(1,1) = v; %
    du(2,1) = (-1/r)*v  - (uu/(1-uu^2))*(v^2-(n^2)/(r^2)) - uu*(1-uu^2);
    du(3,1) = w; % 
    du(4,1) = z*( -((1 + uu^2)/((1-uu^2)^2))*(v^2 - (n^2)/(r^2)) - 1 + 3*uu^2) + w*(-1/r  - ((2*uu)*v)/(1-uu^2)) ;
    
end


