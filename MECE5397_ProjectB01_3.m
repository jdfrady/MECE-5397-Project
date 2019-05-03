%% This code compares ....
%% User Inputs
clc;                          % Clearing the command line 
T = 10;                       % Simulation Time
nt = 20000;                   % Number of time steps
ny = 10;                      % Number of y steps 
nx = 10;                      % Number of x steps
ax = 0;                       % Value of lower x limit
bx = 2*pi;                    % Value of upper x limit
ay = 0;                       % Value of lower y limit
by = 2*pi;                    % Value of upper y limit

%% Definitions and Initializatons
dx = (bx-ax)/(nx+1);          % X step size
dy = (by-ay)/(ny+1);          % Y step Size
dt = T/(nt+1);                % Time step size
u = zeros(nx+2, ny+3, nt+1);  % Allocating memory 
%% Explicit Calculation
% X-axis Boundary Condidions; u(x,y=ay)
temp1 = fb(ay, by);           % Calling function fb for inputs ay & by
temp2 = bx - ax;              % Temporary varibale to simplify for loop
temp3 = gb(ay, by);           % Calling function gb for inputs ay & by
temp4 = temp3 - temp1;        % Temporary varibale to simplify for loop
for t=1:nt+1
    for k = 2:nx+2            % Setting boundary condition 
        u(1,k,t) = temp1 + (((k-1)*dx)/temp2)*temp4; 
    end 
end
% Y-axis Boundary Condidions; u(x=bx,y)
for t=1:nt+1
    for j=2:ny+1              % Setting boundary condition
        u(j,nx+3,t) = gb(((j-1)*dy),by); 
    end
end
% Y-axis Boundary Conditions; u(x=ax,y)
for t=1:nt+1
    for j = 2:ny+1              % Setting boundary condidtion  
        u(j,1,t) = fb(((j-1)*dy),by);
    end
end
% Explicit Calculation
for t = 1:nt                  % Iterating time
    for k = 2:nx              % Iterating in the X direction
        for j = 2:ny+1        % Iterating in the Y direction
            u(j,k,t+1) =  u(j,k,t)+(dt/dx^2)*(u(j,k-1,t)-2*u(j,k,t)+ ...
                u(j,k+1,t))+(dt/dy^2)*(u(j-1,k,t)-2*u(j,k,t)+u(j+1,k,t)); 
        end
    end
    for k = 2:nx+1             % Setting gradient BC (du/dy)|y=by := 0
        u(ny+2,k,t) = u(ny+1,k,t);
    end
end


%% Exact Solution
%% Visualization
    
    surf(u((2:(ny+1)),(2:(nx+1)),nt+1))



