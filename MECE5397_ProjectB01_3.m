%% Main Code
% This code compares the explicit and the implicit approaches to
% demonstrate the accuracy of each method.
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
%% Explicit Calculation
u_explicit = Explicit(T,nt,ny,nx,ax,bx,ay,by);     % Calling the explicit solution
%% Implicit Solution
u_implicit = Implicit_ADI(T,nt,ny,nx,ax,bx,ay,by); % Calling the implict solution
%% Exact Solution
%% Error Calculation

%% Visualization
    figure(1)                                 % Explicit Individual Plot 
    surf(u_explicit((2:(ny+1)),(2:(nx+1)),nt+1))
    
    figure(2)                                 % Implicit Individual Plot
    surf(u_implicit((2:(ny+1)),(2:(nx+1)),nt+1))



