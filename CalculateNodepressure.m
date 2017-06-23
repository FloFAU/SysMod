function [ p_k ] = CalculateNodepressure( deltaP, p1, number_of_knots, knotenindex )
% CalculateNodepressure calculates the pressure at every node.
% The pressure at node 1 (p1) is given (atmosphere)

% INPUT
%   - deltaP (pressure difference)
%   - p1 (presure at node 1) [Pa]
%   - number_of_knots
%   - knotenindex

% OUTPUT
%   - pressure at every node [Pa]


p_k = zeros(number_of_knots,1);

p_k (1) = p1;

    for i = 2:number_of_knots
        for k = 1:number_of_knots 
            if knotenindex((i-1), k) == -1
   
            p_k(i) =abs( p_k(k)-deltaP(i-1)); % Vorsicht bei der Nummerierung der Ströme und Knoten -> Fehler wenn Strom 3 aus Knoten 4 fließen würde
            
            end
        end  
    end
end

