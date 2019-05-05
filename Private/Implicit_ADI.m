%% Implicit Calculation using ADI
%% User Inputs
function u = Implicit_ADI(T,nt,ny,nx,ax,bx,ay,by)
%% Definitions and Initializatons
dx = (bx-ax)/(nx+1);          % X step size
dy = (by-ay)/(ny+1);          % Y step Size
dt = T/(nt+1);                % Time step size
u = zeros(nx+2, ny+2, nt+1);  % Allocating memory 
f = zeros(nx+2,ny+2,nt+1);    % Allocating memory
alpha = zeros(nx+1);          % Allocating memory
g = zeros(nx+1);              % Allocating memory
%% X-axis Boundary Condidions; u(x,y=ay)
temp1 = fb(ay, by);           % Calling function fb for inputs ay & by
temp2 = bx - ax;              % Temporary varibale to simplify for loop
temp3 = gb(ay, by);           % Calling function gb for inputs ay & by
temp4 = temp3 - temp1;        % Temporary varibale to simplify for loop
for t=1:nt+1
    for k = 2:nx            % Setting boundary condition 
        u(1,k,t) = temp1 + (((k-1)*dx)/temp2)*temp4; 
    end 
end
%% Y-axis Boundary Condidions; u(x=bx,y)
for t=1:nt+1
    for j=2:ny+1              % Setting boundary condition
        u(j,nx+2,t) = gb(((j-1)*dy),by); 
    end
end
%% Y-axis Boundary Conditions; u(x=ax,y)
for t=1:nt+1
    for j = 2:ny+1            % Setting boundary condidtion  
        u(j,1,t) = fb(((j-1)*dy),by);
    end
end
%% ADI Method
% To simply the half-step process, I first solved each column at time=t1,
% based on information from time=t0. The second Thomas algorithm solves 
% each row at time=t1 using data from time=t1. This eliminates having to
% create a half step.
lambda1=(dt/(2*dx^2));  % Defining Lambda for the first half-step
lambda2=(dt/(2*dy^2));  % Defining Lambda for the second half-step
a1=(1+2*lambda1);       % Defining matrix variable 'a' for first half step
b1=-lambda1;            % Defining matrix variable 'b' for first half step
c1=-lambda1;            % Defining matrix variable 'c' for first half step
a2=(1+2*lambda2);       % Defining matrix variable 'a' for second half step
b2=-lambda2;            % Defining matrix variable 'b' for second half step
c2=-lambda2;            % Defining matrix variable 'c' for second half step
for t=2:nt+1            % Iterating Time   
    for k=2:ny+1        % First Thomas Algorithm, solving for each row
        for j=2:nx+1
            f(k,j,t)= u(k,j,t-1)+lambda2* ... % Temp variable for half step
                (u(k-1,j,t-1)-2*u(k,j,t-1)+u(k+1,j,t-1));
        end
        alpha(2)=a1;    % This process is manipulating the matrix
        g(2)= f(k,2,t);
        for j=3:nx+1
            alpha(j)=alpha(j)-(b1*c1/alpha(j-1));
            g(j)=f(k,j,t)-((b1*g(j-1))/(alpha(j-1)));
        end
        u(k,nx+1,t)=g(nx+1)/alpha(nx+1);
        for j=1:nx-1   % This is writing the half-steps as whole steps 
            u(k,nx+1-j,t)= (g(nx+1-j)-c1*u(k,nx+2-j,t))/alpha(nx+1-j);
        end
    end
    % Second Thomas Algorithm
    for j=2:nx+1          % Second Half-step
        for k=2:ny+1
            f(k,j,t)= u(k,j,t)+lambda1* ...
                (u(k,j-1,t)-2*u(k,j,t)+u(k,j+1,t));
        end
        alpha(2)=a2;    % Defining alpha(1) for second half step
        g(2)= f(2,j,t);
        for k=3:ny+1
            alpha(k)=alpha(k)-(b2*c2/alpha(k-1));
            g(k)= f(k,j,t)-((b2*g(k-1))/(alpha(k-1)));
        end
        u(ny+1,j,t)=g(ny+1)/alpha(ny+1);
        for k=1:ny-1    % This is overwriting the time step from the 
                        % previous solver
            u(ny+1-k,j,t)= (g(ny+1-k)-c1*u(ny+2-k,j,t))/alpha(ny+1-k);
        end
    end
end
end



