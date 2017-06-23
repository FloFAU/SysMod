%% DESCRIPTION OF THE SYSTEM:

% This approach is based on the method used in "Bestimmung der
% Kühlgasverteilung in Turbogeneratoren durch Kombination der
% Finite-Elemente-Berechnung und der Netzwerkanalyse", Chapter 5 (Zhao)

% All Matrices are simular to the ones used in Zhao's Dissertation, it is
% highly recommended to have look at the Dissertation while working on
% this code.

% The network used for this system is shown on page 69 of the disseration.

% This file is the master-file where the structure of the model is
% determinded.
% To keep the model well structured and easy to work with certain functions
% have been saved in their own files. These function-files have to be put
% in the same folder as this master-file.
% External functional used with this master-file are:

%   - CalculateTemperature.m
%   - CalculateNodepressure.m
%   - CalculateViscosity.m
%   - CalculateDensity.m
%   - CalculateResistance.m


% DATA-INPUT FOR THE SYSTEM:

% There are some informations that have to be typed in the programm each
% time you want to simulate a new network.
% These informations are:

%   - Number of branches, nodes and loops
%   - knotenindex (matrix that contains the connection of the nodes to each other)
%   - B_A (matrix that contains the architecture of the network)
%   - Temperature and pressure of the first node of the system
%   - Length and diameter of the pipes


% (SI-)UNITS USED IN THE PROGRAMM:

% pressure p in Pa
% Temperature T in Kelvin
% length L in m
% diameter d in m
% density rho in kg/m^3
% kinematic viscosity nu in m^2/s
% volume-vlow in ?? (TBD)
% mass-flow in ?? (TBD)

%% This values depend on the system, must be changed every time

b = 6;      % branches
n = 4;      % nodes
c = 3;      % loops (Maschen)

knotenindex = [-1 1 0 0; 0 -1 1 0; 0 -1 0 1 ] ;   % column: node, row: volume-flow
B_A = [-1,0,-1;0,1,-1;1,1,0];                     % loop-matrix for the "real" branches (Äste)

T1=275;          % temperature at node 1
p1=1*10^5;       % pressure at node 1
L = [10;10;10];  % length of tubes 
d = [2;2;2];     % diameter of tubes 
p = zeros(b,1);  
p(1) = 400000;   % sources of the branches

%% Initial values

p_mittel = [300000;290000;290000];                  % has to be automated in the future (TBD)
T_mittel = [222;222;222];                           % has to be automated in the future (TBD)
rho_x = CalculateDensity(p_mittel,T_mittel);        % density
x_old = [12;1;1;133;1;1;1;1;20;2;2;2;1;1;1;1;1;1];  % initial vector, we need to find proper starting values (TBD)

%% Creation of important matrices and vectors

B = zeros(b,c); 
B_S = eye(size(B_A));           % loop-matrix for the "imaginary" branches (Sehnen)
B = [B_A,B_S];                  % combined loop-matrix for all braches (Maschenmatrix)
v_old = x_old((b+1):(2*b));


%% main loop to calculate the solutions-vector of the system of equations

while true
%% Calculating the resistance
nu = CalculateViscosity( p_mittel,T_mittel,rho_x);
Z = CalculateResistance(b,n,c,B_A,rho_x,nu,v_old,L,d);

%% Creating a system of equations
Z_M = B*Z*B';
E_Z = eye(b);
E_K = eye(n-1); % size(E_K)
V_U = zeros(b);
Z_V = zeros(b);
null_one = zeros(b);
null_two = zeros(c);
null_three = zeros(c,b);
middle_A = [E_K,-B_A';null_two,Z_M];
side_A = [null_three;B];
A = [E_Z,-Z,-E_Z;null_one,middle_A,side_A;-V_U,-Z_V,E_Z]; % complete system-matrix

%% Creating the solutions-vector
v_rhs_middle = zeros(n-1,1) ; 
rhs_middle = [v_rhs_middle;B*p] ;
rhs_end = zeros(b,1);
r_k = [-p;rhs_middle;rhs_end];  % complete solutinos-vector

%% Calculation of the new values 
x = linsolve(A,r_k);    % system-solver
x_new = x; 

%% Calculating the pressure at all nodes
kn= [x(1:(n-1))];
p_k = CalculateNodepressure(kn,p1,n,knotenindex);

%% Calculating the average pressure in all pipes:
p_mittel= zeros(n-1,1);

for i = 1:(n-1)
p_mittel(i)= (p_k(i)+p_k(i+1))/2; % wenn wir Druckverluste durch Bauteile in p_Q schreiben, dann: (p_k(1)+p_k(i+1)+p_Q(1)/2)
end

for i = 1:(n-1)
     for k = 1:n 
          if knotenindex((i), k) == -1
           
            p_mittel(i) = (p_k(k)+p_k(i+1))/2;  % hier noch falsch, wenn Quelle vorhanden
            
          end
     end  
 end
  
%% Calculating the temperature at all nodes

T_node = CalculateTemperature(T1,n,knotenindex,p_k);

%% Calculating the average Temperature in all tubes
T_mittel= zeros(n-1,1);

for i = 1:(n-1)
     for k = 1:n 
         if knotenindex((i), k) == -1
           
            T_mittel(i) = (T_node(k)+T_node(i+1))/2; % wenn Bauteil vorhanden falsch; Formel nur für reine Rohre gültig
            
         end
     end  
end

%% Calculating the average density in all tubes
rho_x = CalculateDensity(p_mittel,T_mittel);

%% Calculates the volume-flow
m_real= x_old((b+1):(b+n-1));
v_real = m_real./rho_x;
v_imag = x_old((b+n):(2*b));
v_old = [v_real;v_imag];   
    
%% Check if values are within the tolerance
if  abs(x_old - x_new) <= 0.000001 
   break;
end
    x_old = x_new;
    v_old = x_new((b+1):(2*b));
end

%% Show results:
x_new
Druckdifferenzen = x_new(1:b)
Volumenstroeme = x_new(b+1:2*b)
GesteuerteQuellen = x_new(2*b+1:3*b)