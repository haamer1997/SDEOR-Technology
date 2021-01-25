
global RR;

%Constants
RR = 8.3144621;      % in Pa m3 / mol-K

%for pressure unit converstion
psi2Pa = 6894.757293178;


TinK = 323.15;        % K. Equivalent to 50 Celsius
% TinK = 388.15;   %Bakken oil T

% Number of components
%Data from Chemical Properties Handbook by Yaws, C.L. (1999)
%C1, C2, C3, n-C4, n-C5, n-C6, n-C7, n-C8, n-C9, n-C10, H20, CO2, N2, H2S
pC = [4.604e6 4.88e6 4.249e6 3.797e6 3.369e6 3.012e6 2.736e6 2.486e6 ...
      2.306e6 2.123e6 22.055e6 7.382e6 3.394e6 8.963e6];           % Pa
Tc = [190.58 305.42 369.82 425.18 469.65 507.43 540.26 568.83 595.65 618.45 ...
      647.13 304.19 126.1 373.53];          % deg K
accFact = [1.08e-2 0.099 0.152 0.199 0.249 0.305 0.351 0.396 0.438 ...
           0.484 0.345 0.228 0.040 0.083]; 
MW = [16.043 30.07 44.096 58.123 72.15 86.177 100.204 114.231 128.258 ...
      142.285 18.015 44.01 28.013 34.082]; % in g/mol
Vc = [99.3e-6 147.9e-6 202.9e-6 254.9e-6 312.3e-6 369.9e-6 431.9e-6 ...
      492.1e-6 547.7e-6 603.1e-6 56.0e-6 94.0e-6 90.1e-6 98.5e-6 930.0e-6]; %m^3/mol

Pci = [pC(1) pC(2) pC(10)];
Tci = [Tc(1) Tc(2) Tc(10)];
accFact_i = [accFact(1) accFact(2) accFact(10)];
MWi = [MW(1) MW(2) MW(10)];  %in g/mol
Vci = [Vc(1) Vc(2) Vc(10)];

z1 = 0.85;             % z1 is methane mole fraction
z2 = 0.10;             % z2 is ethane mole fraction
z3 = 1.0 - z1 - z2;
zi = [z1 z2 z3];          % NB z3 is not included as part of the primary variables  

n_c = length(zi);
n_cMinus1 = n_c-1;

deltaij = zeros(3,3);  %It is 3x3 for ternary mixture
deltaij(1,2) = 0.005;
deltaij(2,1) = deltaij(1,2);
deltaij(1,3) = 0.045;
deltaij(3,1) = deltaij(1,3);
deltaij(2,3) = 0.02;
deltaij(3,2) = deltaij(2,3);

%For a gas mixture Ex 4.2.2 on Pg 82 of Multicomponent Mass Trasnfer by Taylor et al.(1993)
D12 = 15.99; %mm^2/s
D13 = 14.43;
D23 = 38.73;

VofC = 15.9;
VofH = 2.31;
VofO = 6.11;
VofN = 4.54;
VofCO2 = 26.7;
VofH20 = 13.1;
VofN2 = 18.5;
VofAir = 19.7;
VofO2 = 16.3;
VofCO = 18.0;

VofC1 = VofC + 4.0*VofH;
VofC2 = 2.0*VofC + 6.0*VofH;
VofC10 = 10.0*VofC + 22.0*VofH;

D12xP = computeGasDijxP(TinK,MWi(1),MWi(2),VofC1,VofC2);
D13xP = computeGasDijxP(TinK,MWi(1),MWi(3),VofC1,VofC10);
D23xP = computeGasDijxP(TinK,MWi(2),MWi(3),VofC2,VofC10);

%specify pressure
pG_inpsi = 4000.0;   %in psia
pG = pG_inpsi*psi2Pa;        %300psi %This is p for the gas phase.

% pG = 47.23e6;    %Initial pore pressure in Pa

D12 = D12xP/pG;
D13 = D13xP/pG;
D23 = D23xP/pG;



B(1,1) = zi(1)/D13 + zi(2)/D12 + zi(3)/D13;
B(1,2) = -zi(1)*(1.0/D12 - 1.0/D13);
B(2,1) = -zi(2)*(1.0/D12 - 1.0/D23);
B(2,2) = zi(1)/D12 + zi(2)/D23 + zi(3)/D23;

Dideal = pinv(B)  %in m^2/s

%Try computing given Dij in Eg 4.2.1
TinK = 273.0;    %in K
MWi = [31.999 28.013 28.01];

D12xP = computeGasDijxP(TinK,MWi(1),MWi(2),VofO2,VofN2);
D13xP = computeGasDijxP(TinK,MWi(1),MWi(3),VofO2,VofCO);
D23xP = computeGasDijxP(TinK,MWi(2),MWi(3),VofN2,VofCO);

%specify pressure
pG = 101.0e3;    %in Pa


D12 = D12xP/pG;
D13 = D13xP/pG;
D23 = D23xP/pG;



B(1,1) = zi(1)/D13 + zi(2)/D12 + zi(3)/D13;
B(1,2) = -zi(1)*(1.0/D12 - 1.0/D13);
B(2,1) = -zi(2)*(1.0/D12 - 1.0/D23);
B(2,2) = zi(1)/D12 + zi(2)/D23 + zi(3)/D23;

Dideal1 = pinv(B); %in m^2/s
Dideal1

% %specify another pressure
% pG_inpsi = 5000.0;   %in psia
% pG = pG_inpsi*psi2Pa;        %300psi %This is p for the gas phase.
% 
% 
% D12 = D12xP/pG;
% D13 = D13xP/pG;
% D23 = D23xP/pG;
% 
% 
% 
% B(1,1) = zi(1)/D13 + zi(2)/D12 + zi(3)/D13;
% B(1,2) = -zi(1)*(1.0/D12 - 1.0/D13);
% B(2,1) = -zi(2)*(1.0/D12 - 1.0/D23);
% B(2,2) = zi(1)/D12 + zi(2)/D23 + zi(3)/D23;
% 
% Dideal2 = pinv(B)  %in m^2/s
% 
% Dideal ./ Dideal2

%computing Fij = I +xi . dln(phi_i)_dxj
% del  = 1e-8;
% 
% [zG, logPhiG] = zFact_n_logPhi(pG,TinK,Pci,Tci, ...
%                          accFact_i, zi, deltaij);
% fGi = zi.*pG.*exp(logPhiG);
% 
% %------------------------------------------------------------------------
% %computing pure gas fugacities...
% [zG0(1), fug0(1)] = zFactPure_n_Fug(pG,TinK,Pci(1),Tci(1),accFact_i(1));
% [zG0(2), fug0(2)] = zFactPure_n_Fug(pG,TinK,Pci(2),Tci(2),accFact_i(2));
% [zG0(3), fug0(3)] = zFactPure_n_Fug(pG,TinK,Pci(3),Tci(3),accFact_i(3));
% 
% fug0 = fug0 .* zi;
% fGi;



%% ------------------------------------------------------------------------

dfGi_dzi = zeros(n_cMinus1,n_cMinus1);
dlnfGi_dzi = zeros(n_cMinus1,n_cMinus1);
dzG_dzi = zeros(1,n_cMinus1);
%little trick to ensure that xi always sums up to 1
zn_decrement = zi;
zn_decrement(n_c) = zi(n_c)-del;
for i = 1:n_cMinus1
%-----------------------
% Single-Phase Gas
%-----------------------
    %resetting 
    zi_incr = zn_decrement;
    %increment/perturb the mole fraction of component i
    zi_incr(i) = zn_decrement(i) + del;
    [zG_atyplusdel, logPhiG_atyplusdel] = zFact_n_logPhi(pG,TinK,Pci,Tci, ...
                         accFact_i, zi_incr, deltaij);
    fGi_atyplusdel_temp = zi_incr.*pG.*exp(logPhiG_atyplusdel);
    fGi_atyplusdel = fGi_atyplusdel_temp(1:n_cMinus1); 
    dfGi_dzi(:,i) = ((fGi_atyplusdel - fGi(1:n_cMinus1)) ./ del)' ;
    dlnfGi_dzi(:,i) = ((log(fGi_atyplusdel) - log(fGi(1:n_cMinus1))) ./ del)' ;

    dzG_dzi(i) = (zG_atyplusdel - zG)/del;
end

%Extract (nc-1) x (nc-1) matrix from dlnfGi_dzi
Fij = dlnfGi_dzi(1:n_c-1,1:n_c-1)
%This is Diffusion coefficient for a real gas
D = Dideal*Fij


% 
% 


  