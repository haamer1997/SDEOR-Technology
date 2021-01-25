function D12xP = computeGasDijxP(TinK,MWi,MWj,Vi,Vj )
%COMPUTEGASDIJXP Summary of this function goes here
%   Detailed explanation goes here
    %Mwi,Mwj in g/mol
    D12xP = 1.013e-2*(TinK^1.75)* ...
        sqrt((MWi+MWj)/(MWi*MWj))/((Vi^(1.0/3.0)+Vj^(1.0/3.0))^2);

end

