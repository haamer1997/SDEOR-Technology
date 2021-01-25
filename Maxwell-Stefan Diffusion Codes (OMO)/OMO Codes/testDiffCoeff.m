


zi = [0.600000	0.399999	0.000001];
ncMinus1 = length(zi)-1;

pG = 2000*6894.757293178;  %4,000 psi

TinK = 323.15;        % K. Equivalent to 50 Celsius
RR = 8.3144621;       % in Pa m3 / mol-K

%Data from Chemical Properties Handbook by Yaws, C.L. (1999)
%C1, C2, C3, n-C4, n-C5, n-C6, n-C7, n-C8, n-C9, n-C10, H20, CO2, N2, H2S
pC = [4.604e6 4.88e6 4.249e6 3.797e6 3.369e6 3.012e6 2.736e6 2.486e6 ...
      2.306e6 2.123e6 22.055e6 7.382e6 3.394e6 8.963e6];           % Pa
Tc = [190.58 305.42 369.82 425.18 469.65 507.43 540.26 568.83 595.65 618.45 ...
      647.13 304.19 126.1 373.53];          % deg K
accFact = [1.08e-2 0.099 0.152 0.199 0.249 0.305 0.351 0.396 0.438 ...
           0.484 0.345 0.228 0.040 0.083]; 
MW = [16.043 30.07 44.096 58.123 72.15 86.177 100.204 114.231 128.258 ...
      142.285 18.015 44.01 28.013 34.082];
Vc = [99.3e-6 147.9e-6 202.9e-6 254.9e-6 312.3e-6 369.9e-6 431.9e-6 ...
      492.1e-6 547.7e-6 603.1e-6 56.0e-6 94.0e-6 90.1e-6 98.5e-6 930.0e-6]; %m^3/mol

Pci = [pC(1) pC(2) pC(3)];
Tci = [Tc(1) Tc(2) Tc(3)];
accFact_i = [accFact(1) accFact(2) accFact(3)];
MWi = [MW(1) MW(2) MW(3)];   
Vci = [Vc(1) Vc(2) Vc(3)];
%         Parachor = [77.0 108.0 150.3];  %Parachors are dimensionless constants

deltaij = zeros(3,3);  %It is 2x2 for binary mixture
deltaij(1,2) = 0.005;
deltaij(2,1) = deltaij(1,2);
deltaij(1,3) = 0.01;
deltaij(3,1) = deltaij(1,3);
deltaij(2,3) = 0.005;
deltaij(3,2) = deltaij(2,3);

VofC = 15.9;
VofH = 2.31;
VofC1 = VofC + 4.0*VofH;
VofC2 = 2.0*VofC + 6.0*VofH;
VofC3 = 3.0*VofC + 8.0*VofH;

D12xP = computeGasDijxP(TinK,MWi(1),MWi(2),VofC1,VofC2);
D13xP = computeGasDijxP(TinK,MWi(1),MWi(3),VofC1,VofC3);
D23xP = computeGasDijxP(TinK,MWi(2),MWi(3),VofC2,VofC3);

Dcoef = computeD(D12xP,D13xP,D23xP,pG,ncMinus1,TinK,Pci,Tci, ...
        accFact_i, zi, deltaij,RR)
    
    
    
    