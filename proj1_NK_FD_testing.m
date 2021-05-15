hold off
clf
clc
clear
hold on


file = fopen("inf.txt", 'w');
for n = [1,2,3] %Here we loop over our hyperparameters and test for
    fprintf(file,'%i\n', n);
    for fn = [1,2,3,4]
    fprintf(file,'%i\n', fn);
    clearvars -except n file fn
    for b = [4:2:16]
%stability i.e. N, b, a,1000 etc.
    a = 0.1; %This is our starting position
    %b = 5; %infity
    N = 50;
    h = (b-a)/N;
    h2 = h*h;
    r = (a:h:b)';
    %u = r/b;
    %u = 1./(1+exp(-r));
    %u = log(r + 1)/log(b+1) -a; 
    %u = -exp(-(r)) + 1 ;
    %u = 1./(1+exp(-r + b/2)); 
    if(fn == 1)
        u = r/b;
    end
    if(fn == 2)
        u = 1./(1+exp(-r + b/2));
    end
    if(fn == 3)
        u = log(r + 1)/log(b+1) -a; 
    end
    if(fn == 4)
        u = -exp(-(r)) + 1;
    end
    rhs=[u(1);u(1:N-1)-2*u(2:N)+u(3:N+1) + h2.*f(r,u,h, N, n);u(N+1)-1];
    tol = 1;
    count = 1;
    while count < 100 && tol > 1e-7
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
        count = count+1; 
    end
        bool = ~isnan(tol) && tol < 1e-7;%tolerence to be considered stable
        fprintf(file,'%f %f\n', b, bool);
    end
        
    end
end %End of for loop
fprintf(file, "TESTING B OVER");
fclose(file);

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