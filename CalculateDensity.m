function [rho_x] = CalculateDensity(pressure,temperature)
% CalculateDensity calculates the density of all pipes
% For the calculations a formular from a scientific paper has been used.

% INPUT:
%   - pressure [Pa]
%   - temperature [K]

% OUTPUT:
%   - density [kg/m^3]


M = 28.949;                       % molar mass [g/mol]
R = 8.3143;                       % molare/universal gas constant [J/mol*K]
rho_x = zeros(size(pressure));

% translation of the units (the following equations only work with Celsius & Mpa)
p = pressure/(10^6);              % Pa -> Mpa
t = temperature-273.15;           % K -> °C

for i=1:size(p)

    % calculates Z with P[MPa] and T[°C]
    Z(i) = 1.00001 - 5.8057*10^(-3)*p(i) + 2.6402*10^(-4)*p(i)^2 - 3.3297*10^(-7)*t(i) + 1.2420*10^(-4)*p(i)*t(i) - 2.0158*10^(-6)*p(i)^2*t(i) + 2.4925*10^(-9)*t(i)^2 - 6.2873*10^(-7)*p(i)*t(i)^2 + 5.4174*10^(-9)*p(i)^2*t(i)^2;

    % calculates density in g/m^3
    rho_x(i) = ((p(i)*10^6*M)/(R*(t(i)+273.15)*Z(i)));

    % calculates density in kg/m^3
    rho_x(i) = rho_x(i)/1000; 
    
end
end

