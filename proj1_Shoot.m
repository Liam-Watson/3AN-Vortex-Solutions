hold off
clf
clc
clear
hold on
%fig = figure("Name", "Nonlinear shooting for a = 0.1, b = 15, tol = 1e-7")
%figure(fig)
a = 0.1; % Boundary one
b = 15; % infinity
n = 2;
%p = 0.2; %Arbitray initial shooting parameter

%file = fopen("inf.txt", 'w');
%fprintf(file,'%6.2f %12.8f\n', n);
for n = [1,2,3]
    ph = 4; %bisection
    pl = -4;
    p = (ph + pl)/2; %bisection
    tol = 1;
    i = 1;
    u0 = [0,p,0,1];
while i < 100 && tol > 1e-7
   [r,u] = ode45(@(r,u) shoot(r,u,n), [a:0.2:b], u0,n);
   %disp(u(end,3))
   %plot(r,u(:,1));
   %pause(0.1)
   %disp((u(end,4)))
   %tol = abs(1-u(end,1));
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
   %disp(p)
   %p = p - ((u(end,1)-1)/(u(end,3)));
   %u0 = [0,p,0,1];
   %disp([p,ph,pl])
   %disp(u0);
   i = i + 1;
   %plot(r,u(:,1))
   %pause(0.1);
end
%    bool = tol < 1e-5;%tolerence to be considered stable
%    fprintf(file,'%6.2f %12.8f\n', b, bool);

plot(r,u(:,1))
legend("n=1", "n=2", "n=3")
title("Nonlinear shooting for a = 0.1, b = 15, tol = 1e-7")
xlabel('r')
ylabel('u')
end

%fclose(file);
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


