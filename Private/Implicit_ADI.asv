%% Implicit Calculation using ADI
%% User Inputs
clc;                          % Clearing the command line 
T = 1;                       % Simulation Time
nt = 2;                   % Number of time steps
ny = 4;                      % Number of y steps 
nx = 4;                      % Number of x steps
ax = 0;                       % Value of lower x limit
bx = 2*pi;                    % Value of upper x limit
ay = 0;                       % Value of lower y limit
by = 2*pi;                    % Value of upper y limit

%% Definitions and Initializatons
dx = (bx-ax)/(nx+1);          % X step size
dy = (by-ay)/(ny+1);          % Y step Size
dt = T/(nt+1);                % Time step size
u = zeros(nx+2, ny+3, nt+1);  % Allocating memory 
f = zeros(nx+1,ny+1,nt+1);    % Allocating memory
alpha = zeros(nx);   %% RECHECK THIS!!!!
g = zeros(nx);
%% X-axis Boundary Condidions; u(x,y=ay)
temp1 = fb(ay, by);           % Calling function fb for inputs ay & by
temp2 = bx - ax;              % Temporary varibale to simplify for loop
temp3 = gb(ay, by);           % Calling function gb for inputs ay & by
temp4 = temp3 - temp1;        % Temporary varibale to simplify for loop
for t=1:nt+1
    for k = 2:nx+2            % Setting boundary condition 
        u(1,k,t) = temp1 + (((k-1)*dx)/temp2)*temp4; 
    end 
end
%% Y-axis Boundary Condidions; u(x=bx,y)
for t=1:nt+1
    for j=2:ny+1              % Setting boundary condition
        u(j,nx+3,t) = gb(((j-1)*dy),by); 
    end
end

%% Y-axis Boundary Conditions; u(x=ax,y)
for t=1:nt+1
    for j = 2:ny+1            % Setting boundary condidtion  
        u(j,1,t) = fb(((j-1)*dy),by)
    end
end

%% First half-step
lambda1=dt/(2*dx^2);
lambda2=dt/(2*dy^2);
a=-lambda1;
b=(1+2*lambda1);
c=-lambda1;
for t=1:1
    for k=2:ny
        for j=2:nx
            f(k,j,t)= u(k,j,t)+lambda2*(u(k-1,j,t)-2*u(k,j,t)+u(k+1,j,t));
        end
    end
    % Thomas Algorithm
    for k=2:ny
        alpha(1)=a;
        g(1)= f(k,2,t);
        for j=2:nx
            alpha(j)=alpha(j)-(b*c/alpha(j-1));
            g(j)=f(k,j,t)-((b*g(j-1))/(alpha(j-1)));
        end
        un=g(nx)/alpha(nx);
        for j=1:nx-1
            u(k,nx-j,t)= (g(nx-j)-c*u(k,nx-j+1,t))/alpha(nx-j)
        end
    end
    
    
    
    
    
    
    
    
    
end






