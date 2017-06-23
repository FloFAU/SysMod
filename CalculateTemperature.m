function [T_node] = CalculateTemperature(T1,number_of_nodes,knotenindex,p_node)
% CalculateTemperature calculates the temperature at every node (in Kelvin)
% assumption: change of the different states happens reversibel adiabat

% INPUT
%   - pemperature at first node T1 [K]
%   - number_of_nodes 
%   - knotenindex
%   - pressure at the nodes p_node [Pa]

% OUTPUT
%   - temperatrure at the nodes [K]


T_node = zeros(number_of_nodes,1);
T_node (1) = T1; % in K

cp = 1.005;  % specific heat capacity (isobar) for air [kJ/(kg*K)] 
cv = 0.718;  % specific heat capacity (isochor) for air [kJ/(kg*K)] 

n = cp/cv; % assumption: reversibel adiabat

    for i = 2:number_of_nodes
        for k = 1:number_of_nodes 
            if knotenindex((i-1), k) == -1
   
                T_node(i) = (T_node(k)*(p_node(i)/p_node(k))^((n-1)/n)); % Vorsicht bei der Nummerierung der Ströme und Knoten -> Fehler wenn Strom 3 aus Knoten 4 fließen würde
            
            end
        end  
    end
end

