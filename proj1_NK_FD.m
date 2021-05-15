
hold off
clf
clc
clear
hold on

%h = (b-a)/N;

%r = (a:h:b)';

%Initial guesses
%u = r/b;


%n = 1; %vorticity 
%file = fopen("inf.txt", 'w');
for n = [1,2,3] %Here we loop over our hyperparameters and test for
%stability i.e. N, b, a,1000 etc.
    a = 0; %This is our starting position
    b = 7; %infity
    N = 50;
    if n == 1
        a = 0.1;
    end
    h = (b-a)/N;
    h2 = h*h;
    r = (a:h:b)';
    %u = r/b;
    %u = 1./(1+exp(-r));
    u = log(r + 1)/log(b+1) -a; 
    %u = -exp(-(r)) + 1 ;
    
    %u = 1./(1+exp(-r + b/2)); 
    rhs=[u(1);u(1:N-1)-2*u(2:N)+u(3:N+1) + h2.*f(r,u,h, N, n);u(N+1)-1];
    tol = 1;
    count = 1;
    %n = 2
    while count < 100000 && tol > 1e-7
        dd = [0; -2 - h2.*q(r,u,h,N, n) ; 0];
        upper = diag([0;1 - (h*p(r,u,h,N))./2],1);
        lower = diag([1 + (h*p(r,u,h,N))./2;0],-1);
        J = diag(dd,0) + upper + lower;  %Set three diagonal
        J(1,1)=1; %J(1,2)=0; J(1,3)=0; %Set BC 1
        J(N+1,N+1)=1; %J(N+1,N)=0; J(N+1,N-1)=0; %Set BC 2
        z = -J\rhs; %Solve for correction
        u = u + z; 
        rhs=[u(1);(u(1:N-1)-2*u(2:N)+u(3:N+1)) + h2.*f(r,u,h,N,n);u(N+1)-1]; 
        tol = norm(rhs); %Display tolerance
        %plot(r,u)
        %pause(0.05);
        count = count+1; 
    end
%    bool = tol < 1e-5;%tolerence to be considered stable
%    fprintf(file,'%6.2f %12.8f\n', b, bool); %Write stability data
    disp(tol)
    plot(r,u)
    legend("n=1", "n=2", "n=3")
    title("NK for b = 7, tol = 1e-7, logarithm guess, N = 50")
    xlabel('r')
    ylabel('u')
end %End of for loop
%fclose(file);500
%This is f in y'' = f
function [fu] = f(r, u, h, N, n)
    uBefore = u(1:N-1);
    uAfter = u(3:N+1);
    u = u(2:N);
    fu = 1*((1./r(2:N)).*((uAfter - uBefore)/(2*h)) + ((u)./(1-u.^2)).*(((uAfter - uBefore)/(2*h)).^2 - n^2./(r(2:N).^2)) + u.*(1-u.^2));
end 
%This is p in a linear ODE
function [p] = p(r,u,h, N)
    uBefore = u(1:N-1);
    uAfter = u(3:N+1);
    uu = u(2:N);
    r = r(2:N);
    p = -1*(1./r + (2*uu.*((uAfter-uBefore)/(2*h)))./(1-uu.^2));
end


%q(r)
function [q] =  q(r,u, h , N, n)
    r = r(2:N);
    uBefore = u(1:N-1);
    uAfter = u(3:N+1);
    uu = u(2:N);
    q= -1*(-3*uu.^2 + 1 + 1./((1-uu.^2)).*(((uAfter-uBefore)/(2*h)).^2 - (n^2)./(r.^2) )); 
end