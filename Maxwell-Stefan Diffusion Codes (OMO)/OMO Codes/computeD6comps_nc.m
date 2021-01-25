function D = computeD6comps_nc(fluid,state,Constants)
%COMPUTED Summary of this function goes here
%   Detailed explanation goes here
%    
    %relative to total mass of hydrocarbon
    B = zeros(fluid.ncMinus1,fluid.ncMinus1);
    
    %relative to heaviest component
    B(1,1) =  state.zi_init(2)/fluid.D12xP + state.zi_init(3)/fluid.D13xP + state.zi_init(4)/fluid.D14xP + state.zi_init(5)/fluid.D15xP + state.zi_init(6)/fluid.D16xP; % pg 87
    B(1,2) = -state.zi_init(1)/fluid.D12xP;
    B(1,3) = -state.zi_init(1)/fluid.D13xP;
    B(1,4) = -state.zi_init(1)/fluid.D14xP;
    B(1,5) = -state.zi_init(1)/fluid.D15xP;
    B(2,1) = -state.zi_init(2)/fluid.D12xP;
    B(2,2) =  state.zi_init(1)/fluid.D12xP + state.zi_init(3)/fluid.D23xP + state.zi_init(4)/fluid.D24xP + state.zi_init(5)/fluid.D25xP + state.zi_init(6)/fluid.D26xP;
    B(2,3) = -state.zi_init(2)/fluid.D23xP;
    B(2,4) = -state.zi_init(2)/fluid.D24xP;
    B(2,5) = -state.zi_init(2)/fluid.D25xP;
    B(3,1) = -state.zi_init(3)/fluid.D13xP;
    B(3,2) = -state.zi_init(3)/fluid.D23xP;
    B(3,3) =  state.zi_init(1)/fluid.D13xP + state.zi_init(2)/fluid.D23xP + state.zi_init(4)/fluid.D34xP + state.zi_init(5)/fluid.D35xP + state.zi_init(6)/fluid.D36xP;
    B(3,4) = -state.zi_init(3)/fluid.D34xP;
    B(3,5) = -state.zi_init(3)/fluid.D35xP;
    B(4,1) = -state.zi_init(4)/fluid.D14xP;
    B(4,2) = -state.zi_init(4)/fluid.D24xP;
    B(4,3) = -state.zi_init(4)/fluid.D34xP;
    B(4,4) =  state.zi_init(1)/fluid.D14xP + state.zi_init(2)/fluid.D24xP + state.zi_init(3)/fluid.D34xP + state.zi_init(5)/fluid.D45xP + state.zi_init(6)/fluid.D46xP;
    B(4,5) = -state.zi_init(4)/fluid.D45xP;
    B(5,1) = -state.zi_init(5)/fluid.D15xP;
    B(5,2) = -state.zi_init(5)/fluid.D25xP;
    B(5,3) = -state.zi_init(5)/fluid.D35xP;
    B(5,4) = -state.zi_init(5)/fluid.D45xP;
    B(5,5) =  state.zi_init(1)/fluid.D15xP + state.zi_init(2)/fluid.D25xP + state.zi_init(3)/fluid.D35xP + state.zi_init(4)/fluid.D45xP + state.zi_init(6)/fluid.D56xP;
    
    Dideal = pinv(B)./state.p_init*1e6;  %in m^2/s       
        
    del  = 1e-8;
    [~, logPhiG] = zFact_n_logPhi(state.p_init*1e6,state.T,fluid.Pci*1e6,fluid.Tci, ...
                             fluid.accFact_i, state.zi_init, fluid.deltaij,Constants.RR);

    dlnPhi_i_dzj = zeros(fluid.ncMinus1,fluid.ncMinus1);
    %little trick to ensure that xi always sums up to 1
    zn_decrement = state.zi_init;
    zn_decrement(fluid.n_c) = state.zi_init(fluid.n_c)-del;
    for i = 1:fluid.ncMinus1
        %resetting 
        zi_incr = zn_decrement;
        %increment/perturb the mole fraction of component i
        zi_incr(i) = zn_decrement(i) + del;
        [~, logPhiG_atyplusdel] = zFact_n_logPhi(state.p_init*1e6,state.T,fluid.Pci*1e6,fluid.Tci, ...
                             fluid.accFact_i, zi_incr, fluid.deltaij,Constants.RR);
        dlnPhi_i_dzj(:,i) = state.zi_init(i)*((logPhiG_atyplusdel(1:fluid.ncMinus1) - logPhiG(1:fluid.ncMinus1)) ./ del)';
    end

    %Extract (nc-1) x (nc-1) matrix from dlnfGi_dzi
    Fij = eye(fluid.ncMinus1) + dlnPhi_i_dzj; %note that dlnPhi_i_dzj is already multiplied by zi
    %This is Diffusion coefficient for a real gas
    D = Dideal*Fij;

end

