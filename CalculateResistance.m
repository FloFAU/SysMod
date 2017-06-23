function [Z] = CalculateResistance(b,n,c,B_A,rho,nu,v_old,L,d)
% CalculateResistance calculates the resistance of all branches

% INPUT:
%   - branches b
%   - nodes n
%   - loops c
%   - B_A (loop-matrix, Äste)
%   - density rho [kg/m^3]
%   - v_old (inital solutions-vector)
%   - Length of pipe L [m]
%   - diameter of pipe d [m]

% OUTPUT:
%   - resistance matrix Z

B_S = eye(size(B_A));

% Calculating the resistance of all "real" branches (Äste)
R_e = zeros(n-1,1);
z_a = zeros(n-1,1);
lambda = zeros((n-1),1);

for i=1:(n-1) 
    R_e(i) = (4/(pi*d(i)*nu(i)))*v_old(i);  % calcultaion of the Reynols-Number
    
    if R_e(i) < 2320                        % calculation of lambda (3 different equations/cases)
        lambda(i) = 64/R_e(i);    
    elseif 2320<=R_e(i) < 10^5 
        lambda(i) = 0.3164*R_e(i).^(-0.25);    
    else 
        lambda(i) = 0.0032 + (0.221*R_e(i).^(-0.237));
    end
    
    z_a(i) = lambda(i) *(8*rho(i)/pi.^2)*(L(i)/d(i).^5)*v_old(i); 
end

% Resistance values for the "imaginary" branches (Sehnenwiderstände Z_S)
Z_A=zeros(n-1); % To calculate the resistance of the "imaginary" branches (Sehnen) a resistance-Matrix of branches (Äste) is needed
for i=1:(n-1)
        Z_A(i,i) = z_a(i);
end

Zaehler = B_A*Z_A*v_old(1:(n-1)); % fehlen Druckquellen und gesteuerte Druckabfälle: -B_A*P_q(1:(n-1)) +B_A*P_st_old(1:(n-1))
Nenner = B_S*v_old((n):b);
Z_s = zeros(c,1); 
for i=1:c
   Z_s(i)= Zaehler(i)/Nenner(i);  
end

% Adding the resistance of Äste and Sehnen to get a matrix with the resistance of all branches
Z_as = [z_a; Z_s];
Z = zeros(b);

for i = 1:b
    Z(i,i) = Z_as(i);
end
end

