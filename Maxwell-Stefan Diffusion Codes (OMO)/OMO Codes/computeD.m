function D = computeD(D12xP,D13xP,D23xP,pG,n_cMinus1,TinK,Pci,Tci, ...
                             accFact_i, zi, deltaij,RR)
%COMPUTED Summary of this function goes here
%   Detailed explanation goes here
    n_c = n_cMinus1+1;
    
    B = zeros(n_cMinus1,n_cMinus1);
    B(1,1) = zi(1)/D13xP + zi(2)/D12xP + zi(3)/D13xP;
    B(1,2) = -zi(1)*(1.0/D12xP - 1.0/D13xP);
    B(2,1) = -zi(2)*(1.0/D12xP - 1.0/D23xP);
    B(2,2) = zi(1)/D12xP + zi(2)/D23xP + zi(3)/D23xP;
%     B(1,1) = zi(2)/D12xP + zi(3)/D13xP;
%     B(1,2) = -zi(1)*(1.0/D12xP);
%     B(2,1) = -zi(2)*(1.0/D12xP);
%     B(2,2) = zi(1)/D12xP + zi(3)/D23xP;
    
    Dideal = pinv(B)./pG;  %in m^2/s       
        
    del  = 1e-8;
    [~, logPhiG] = zFact_n_logPhi(pG,TinK,Pci,Tci, ...
                             accFact_i, zi, deltaij,RR);

    dlnPhi_i_dzj = zeros(n_cMinus1,n_cMinus1);
    %little trick to ensure that xi always sums up to 1
    zn_decrement = zi;
    zn_decrement(n_c) = zi(n_c)-del;
    for i = 1:n_cMinus1
        %resetting 
        zi_incr = zn_decrement;
        %increment/perturb the mole fraction of component i
        zi_incr(i) = zn_decrement(i) + del;
        [~, logPhiG_atyplusdel] = zFact_n_logPhi(pG,TinK,Pci,Tci, ...
                             accFact_i, zi_incr, deltaij,RR);
        dlnPhi_i_dzj(:,i) = zi(i)*((logPhiG_atyplusdel(1:n_cMinus1) - logPhiG(1:n_cMinus1)) ./ del)';
    end

    %Extract (nc-1) x (nc-1) matrix from dlnfGi_dzi
    Fij = eye(n_cMinus1) + dlnPhi_i_dzj; %note that dlnPhi_i_dzj is already multiplied by zi
    %This is Diffusion coefficient for a real gas
    D = Dideal*Fij;

end

