function Dideal = dIdeal(state,fluidMixture,rel2totalMass)
%gasDijxP Summary of this function goes here
%   Detailed explanation goes here
    %Mwi,Mwj in g/mol
    
    nc = numel(state.components);
    %From Taylor and Krishna pg 69
    VofC = 15.9;
    VofH = 2.31;
    %         VofO = 6.11;
    VofC1 = VofC + 4.0*VofH;
    VofC2 = 2.0*VofC + 6.0*VofH;
    %         VofCO2 = VofC + 2.0*VofO;
    %         VofCO2 = 26.7;
    VofC3 = 3.0*VofC + 8.0*VofH;
    VofC4 = 4.0*VofC + 10.0*VofH;
    VofC5 = 5.0*VofC + 12.0*VofH;
    VofC6 = 6.0*VofC + 14.0*VofH;

    B = zeros(nc-1,nc-1);
    switch nc
        case 2
            D12xP = gasDijxP(state.T,fluidMixture.molarMass(1),fluidMixture.molarMass(2),VofC1,VofC2);
            Dideal = D12xP;  % TO DO: This can't be right. Focus is on nc>2 for now
            
        case 3
            D12xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(3),VofC1,VofC2); %note the alteration of indices bcos C1 is 2nd component
            D13xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(1),VofC1,VofC3);
            D23xP = gasDijxP(state.T,fluidMixture.molarMass(3),fluidMixture.molarMass(1),VofC2,VofC3);
 
            if(rel2totalMass)
                %relative to total mass of hydrocarbon  pg 87
                B(1,1) = state.y(2)/D13xP + state.y(3)/D12xP + state.y(1)/D13xP; %note the alteration of indices bcos C1 is 2nd component in fluid mixture
                B(1,2) = -state.y(2)*(1.0/D12xP - 1.0/D13xP);
                B(2,1) = -state.y(3)*(1.0/D12xP - 1.0/D23xP);
                B(2,2) = state.y(2)/D12xP + state.y(3)/D23xP + state.y(1)/D23xP;
            else
                %relative to heaviest component 
                B(1,1) = state.y(3)/D12xP + state.y(1)/D13xP;
                B(1,2) = -state.y(2)*(1.0/D12xP);
                B(2,1) = -state.y(3)*(1.0/D12xP);
                B(2,2) = state.y(2)/D12xP + state.y(1)/D23xP;
            end
    
            Dideal = pinv(B)./state.pressure;  %in m^2/s 

        case 4
            D12xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(3),VofC1,VofC2);
            D13xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(4),VofC1,VofC3);
            D14xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(1),VofC1,VofC4);
            D23xP = gasDijxP(state.T,fluidMixture.molarMass(3),fluidMixture.molarMass(4),VofC2,VofC3);
            D24xP = gasDijxP(state.T,fluidMixture.molarMass(3),fluidMixture.molarMass(1),VofC2,VofC4);
            D34xP = gasDijxP(state.T,fluidMixture.molarMass(4),fluidMixture.molarMass(1),VofC3,VofC4);
            
            if(rel2totalMass)
                %relative to total mass of hydrocarbon  pg 87
                B(1,1) = state.y(2)/D14xP + state.y(3)/D12xP + state.y(4)/D13xP + state.y(1)/D14xP; % pg 87
                B(1,2) = -state.y(2)*(1.0/D12xP - 1.0/D14xP);
                B(1,3) = -state.y(2)*(1.0/D13xP - 1.0/D14xP);
                B(2,1) = -state.y(3)*(1.0/D12xP - 1.0/D24xP);
                B(2,2) = state.y(2)/D12xP + state.y(3)/D24xP + state.y(4)/D23xP + state.y(1)/D24xP;
                B(2,3) = -state.y(3)*(1.0/D23xP - 1.0/D24xP);
                B(3,1) = -state.y(4)*(1.0/D13xP - 1.0/D34xP);
                B(3,2) = -state.y(4)*(1.0/D23xP - 1.0/D34xP);%todo
                B(3,3) = state.y(2)/D13xP + state.y(3)/D23xP + state.y(4)/D34xP + state.y(1)/D34xP;
            else
                %relative to heaviest component 
                
            end

            Dideal = pinv(B)/state.pressure;  %in m^2/s 

        case 5
            D12xP = gasDijxP(state.T,fluidMixture.molarMass(1),fluidMixture.molarMass(2),VofC1,VofC2);
            D13xP = gasDijxP(state.T,fluidMixture.molarMass(1),fluidMixture.molarMass(3),VofC1,VofC3);
            D14xP = gasDijxP(state.T,fluidMixture.molarMass(1),fluidMixture.molarMass(4),VofC1,VofC4);
            D15xP = gasDijxP(state.T,fluidMixture.molarMass(1),fluidMixture.molarMass(5),VofC1,VofC5);
            D23xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(3),VofC2,VofC3);
            D24xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(4),VofC2,VofC4);
            D25xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(5),VofC2,VofC5);
            D34xP = gasDijxP(state.T,fluidMixture.molarMass(3),fluidMixture.molarMass(4),VofC3,VofC4);
            D35xP = gasDijxP(state.T,fluidMixture.molarMass(3),fluidMixture.molarMass(5),VofC3,VofC5);
            D45xP = gasDijxP(state.T,fluidMixture.molarMass(4),fluidMixture.molarMass(5),VofC4,VofC5);
            
            if(rel2totalMass)
                %relative to total mass of hydrocarbon  pg 87
                B(1,1) = state.y(2)/D15xP + state.y(3)/D12xP + state.y(4)/D13xP + state.y(5)/D14xP + state.y(1)/D15xP; % pg 87
                B(1,2) = -state.y(2)*(1.0/D12xP - 1.0/D15xP);
                B(1,3) = -state.y(2)*(1.0/D13xP - 1.0/D15xP);
                B(1,4) = -state.y(2)*(1.0/D14xP - 1.0/D15xP);
                B(2,1) = -state.y(3)*(1.0/D12xP - 1.0/D25xP);
                B(2,2) = state.y(2)/D12xP + state.y(3)/D25xP + state.y(4)/D23xP + state.y(5)/D24xP + state.y(1)/D25xP;
                B(2,3) = -state.y(3)*(1.0/D23xP - 1.0/D25xP);
                B(2,4) = -state.y(3)*(1.0/D24xP - 1.0/D25xP);
                B(3,1) = -state.y(4)*(1.0/D13xP - 1.0/D35xP);
                B(3,2) = -state.y(4)*(1.0/D23xP - 1.0/D35xP);
                B(3,3) = state.y(2)/D13xP + state.y(3)/D23xP + state.y(4)/D35xP + state.y(5)/D34xP + state.y(1)/D35xP;
                B(3,4) = -state.y(4)*(1.0/D34xP - 1.0/D35xP);
                B(4,1) = -state.y(5)*(1.0/D14xP - 1.0/D45xP);
                B(4,2) = -state.y(5)*(1.0/D24xP - 1.0/D45xP);
                B(4,3) = -state.y(5)*(1.0/D34xP - 1.0/D45xP);
                B(4,4) = state.y(2)/D14xP + state.y(3)/D24xP + state.y(4)/D34xP + state.y(5)/D45xP + state.y(1)/D45xP;  
            else
                %relative to heaviest component 
                B(1,1) =  state.y(3)/D12xP + state.y(4)/D13xP + state.y(5)/D14xP + state.y(1)/D15xP; % pg 87
                B(1,2) = -state.y(1)/D12xP;
                B(1,3) = -state.y(1)/D13xP;
                B(1,4) = -state.y(1)/D14xP;
                B(2,1) = -state.y(2)/D12xP;
                B(2,2) = state.y(2)/D12xP + state.y(4)/D23xP + state.y(5)/D24xP + state.y(1)/D25xP;
                B(2,3) = -state.y(2)/D23xP;
                B(2,4) = -state.y(2)/D24xP;
                B(3,1) = -state.y(3)/D13xP;
                B(3,2) = -state.y(3)/D23xP;
                B(3,3) = state.y(2)/D13xP + state.y(3)/D23xP + state.y(5)/D34xP + state.y(1)/D35xP;
                B(3,4) = -state.y(3)/D34xP;
                B(4,1) = -state.y(4)/D14xP;
                B(4,2) = -state.y(4)/D24xP;
                B(4,3) = -state.y(4)/D34xP;
                B(4,4) = state.y(2)/D14xP + state.y(3)/D24xP + state.y(4)/D34xP + state.y(1)/D45xP;
            end

            Dideal = pinv(B)/state.pressure;  %in m^2/s 
            
        case 6
            D12xP = gasDijxP(state.T,fluidMixture.molarMass(1),fluidMixture.molarMass(2),VofC1,VofC2);
            D13xP = gasDijxP(state.T,fluidMixture.molarMass(1),fluidMixture.molarMass(3),VofC1,VofC3);
            D14xP = gasDijxP(state.T,fluidMixture.molarMass(1),fluidMixture.molarMass(4),VofC1,VofC4);
            D15xP = gasDijxP(state.T,fluidMixture.molarMass(1),fluidMixture.molarMass(5),VofC1,VofC5);
            D16xP = gasDijxP(state.T,fluidMixture.molarMass(1),fluidMixture.molarMass(6),VofC1,VofC6);
            D23xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(3),VofC2,VofC3);
            D24xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(4),VofC2,VofC4);
            D25xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(5),VofC2,VofC5);
            D26xP = gasDijxP(state.T,fluidMixture.molarMass(2),fluidMixture.molarMass(6),VofC2,VofC6);
            D34xP = gasDijxP(state.T,fluidMixture.molarMass(3),fluidMixture.molarMass(4),VofC3,VofC4);
            D35xP = gasDijxP(state.T,fluidMixture.molarMass(3),fluidMixture.molarMass(5),VofC3,VofC5);
            D36xP = gasDijxP(state.T,fluidMixture.molarMass(3),fluidMixture.molarMass(6),VofC3,VofC6);
            D45xP = gasDijxP(state.T,fluidMixture.molarMass(4),fluidMixture.molarMass(5),VofC4,VofC5);
            D46xP = gasDijxP(state.T,fluidMixture.molarMass(4),fluidMixture.molarMass(6),VofC4,VofC6);
            D56xP = gasDijxP(state.T,fluidMixture.molarMass(5),fluidMixture.molarMass(6),VofC5,VofC6);
            
            if(rel2totalMass)
                %relative to total mass of hydrocarbon  pg 87
                
            else
                %relative to heaviest component
                B(1,1) =  state.y(3)/D12xP + state.y(4)/D13xP + state.y(5)/D14xP + state.y(6)/D15xP + state.y(1)/D16xP; % pg 87
                B(1,2) = -state.y(2)/D12xP;
                B(1,3) = -state.y(2)/D13xP;
                B(1,4) = -state.y(2)/D14xP;
                B(1,5) = -state.y(2)/D15xP;
                B(2,1) = -state.y(3)/D12xP;
                B(2,2) =  state.y(2)/D12xP + state.y(4)/D23xP + state.y(5)/D24xP + state.y(6)/D25xP + state.y(1)/D26xP;
                B(2,3) = -state.y(3)/D23xP;
                B(2,4) = -state.y(3)/D24xP;
                B(2,5) = -state.y(3)/D25xP;
                B(3,1) = -state.y(4)/D13xP;
                B(3,2) = -state.y(4)/D23xP;
                B(3,3) =  state.y(2)/D13xP + state.y(3)/D23xP + state.y(5)/D34xP + state.y(6)/D35xP + state.y(1)/D36xP;
                B(3,4) = -state.y(4)/D34xP;
                B(3,5) = -state.y(4)/D35xP;
                B(4,1) = -state.y(5)/D14xP;
                B(4,2) = -state.y(5)/D24xP;
                B(4,3) = -state.y(5)/D34xP;
                B(4,4) =  state.y(2)/D14xP + state.y(3)/D24xP + state.y(4)/D34xP + state.y(6)/D45xP + state.y(1)/D46xP;
                B(4,5) = -state.y(5)/D45xP;
                B(5,1) = -state.y(6)/D15xP;
                B(5,2) = -state.y(6)/D25xP;
                B(5,3) = -state.y(6)/D35xP;
                B(5,4) = -state.y(6)/D45xP;
                B(5,5) =  state.y(2)/D15xP + state.y(3)/D25xP + state.y(4)/D35xP + state.y(5)/D45xP + state.y(1)/D56xP;
            end
    
            Dideal = pinv(B)/state.pressure;  %in m^2/s                    

        otherwise
            disp('Number of components is not between 2 and 6 (inclusive)')
    end

end

function D12xP = gasDijxP(TinK,MWi,MWj,Vi,Vj )
%gasDijxP Summary of this function goes here
%   Detailed explanation goes here
    %Mwi,Mwj in g/mol
    MWi = MWi*1e3;
    MWj = MWj*1e3;
    D12xP = 1.013e-2*(TinK^1.75)* ...
        sqrt((MWi+MWj)/(MWi*MWj))/((Vi^(1.0/3.0)+Vj^(1.0/3.0))^2);
end

