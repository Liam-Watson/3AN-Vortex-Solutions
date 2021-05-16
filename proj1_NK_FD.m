hold off
clf
clc
clear

for n = [1,2,3] % Loop over all n=1,2,3 and obtain numerical solutions
    a = 0; %The left BC
    b = 4; %infity (the right BC)
    N = 50; %Define number of mesh points
    h = (b-a)/N; %Define the mesh step size according to scheme
    h2 = h*h; %Simplification of code
    r = (a:h:b)'; %Set up r mesh with step size h between a and b
    u = r/b; %Initial guess is linear, other guesses are bellow
    %u = log(r + 1)/log(b+1) -a; 
    %u = -exp(-(r)) + 1 ;
    %u = 1./(1+exp(-r + b/2)); 
    rhs=[u(1);u(1:N-1)-2*u(2:N)+u(3:N+1) + h2.*f(r,u,h, N, n);u(N+1)-1]; %Set initial RHS based on guess
    tol = 1; %Initialize tolerance value
    count = 1; %Initialize iteration counter
    while count < 100 && tol > 1e-7
        dd = [0; -2 - h2.*q(r,u,h,N, n) ; 0];
        upper = diag([0;1 - (h*p(r,u,h,N))./2],1);
        lower = diag([1 + (h*p(r,u,h,N))./2;0],-1);
        J = diag(dd,0) + upper + lower;  %Set three diagonal
        J(1,1)=1; %Set BC 1
        J(N+1,N+1)=1; %Set BC 2
        z = -J\rhs; %Solve for correction
        u = u + z; %Update u with new correction
        rhs=[u(1);(u(1:N-1)-2*u(2:N)+u(3:N+1)) + h2.*f(r,u,h,N,n);u(N+1)-1]; %Update RHS values
        tol = norm(rhs); %Display tolerance
        fig1 = figure(1);
        plot(r,u); %Plot updated curve of u
        hold off;
        pause(0.1);
        count = count+1; %Increment iteration counter
    end
    %Plot final solution u for each n
    fig2 = figure(2);
    plot(r,u);
    hold on;
    legend("n=1", "n=2", "n=3");
    title("NK for b = 7, tol = 1e-7, logarithm guess, N = 50");
    xlabel('r');
    ylabel('u');
    pause(1);
end %End of n for loop

%This is f in y'' = f
function [fu] = f(r, u, h, N, n)
    uBefore = u(1:N-1); %This is u_{j-1}
    uAfter = u(3:N+1); %This is u_{j+1}
    u = u(2:N); %This is u_{j}
    %Evaluate and return f(2:N). Note the FD approximations of u' 
    fu = (1./r(2:N)).*((uAfter - uBefore)/(2*h)) + ((u)./(1-u.^2)).*(((uAfter - uBefore)/(2*h)).^2 - n^2./(r(2:N).^2)) + u.*(1-u.^2);
end 

%This is p(r) in a linear kantorovich equations
function [p] = p(r,u,h, N)
    uBefore = u(1:N-1); %This is u_{j-1}
    uAfter = u(3:N+1); %This is u_{j+1}
    uu = u(2:N); %This is u_{j}
    r = r(2:N); %Find r excluding the end points
    %Evaluate and return p(2:N)
    p = -1*(1./r + (2*uu.*((uAfter-uBefore)/(2*h)))./(1-uu.^2));
end

%This is q(r) in the linear kantorovich equations
function [q] =  q(r,u, h , N, n)
    r = r(2:N); %Find r excluding the end points
    uBefore = u(1:N-1);  %This is u_{j-1}
    uAfter = u(3:N+1); %This is u_{j+1}
    uu = u(2:N); %This is u_{j}
    %Evaluate and return q(2:N)
    q= -1*(-3*uu.^2 + 1 + 1./((1-uu.^2)).*(((uAfter-uBefore)/(2*h)).^2 - (n^2)./(r.^2) )); 
end