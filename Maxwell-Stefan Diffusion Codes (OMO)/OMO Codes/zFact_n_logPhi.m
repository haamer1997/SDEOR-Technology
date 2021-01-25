% '------------------------------------------------------------------------
% 'Multicomponent Z Factor using Peng_Robinson Equation of States
% '------------------------------------------------------------------------
function [z, logPhi] = zFact_n_logPhi(pres,TinK,pCrit,tCrit,accFact, x, deltaij,RR)
    Tr     = TinK./tCrit;
    kappa  = (-0.26992.*accFact + 1.54226).*accFact + 0.37464;
    kappa(kappa>1.07829) = ((0.01667.*accFact(accFact>0.5)-0.1644) ...
        .*accFact(accFact>0.5)+1.485).*accFact(accFact>0.5)+0.3796;

    sqrtAlpha = (1.0 + kappa .* (1.0-sqrt(Tr)));
    alpha   = sqrtAlpha .* sqrtAlpha ;
    aa      = alpha .* (0.457235*RR*RR) .* tCrit .* tCrit ./ pCrit;
    bb      = 0.077796 .* RR .* tCrit ./ pCrit;
    
    bMix   = x*bb';
    aij = (1-deltaij).*sqrt(aa'*aa); 
    
    aMix = sum(sum((x'*x).*aij));
    
    A_ofP = aMix * pres / (RR * RR * TinK * TinK);
    B_ofP = bMix * pres/(RR * TinK);

    a2 = -(1.0d0 - B_ofP);
    a1 =  A_ofP - 3.0d0*B_ofP*B_ofP - 2.0d0*B_ofP;
    a0 = -B_ofP*( A_ofP - B_ofP*(1.0d0 + B_ofP) );
    
    [myRoots, numRoots] = rootsOfCubicEqn( a2, a1, a0 );
    
    if (~isreal([myRoots, B_ofP ]))
        myRoots
    end
    z = 0.0;
    if (numRoots==1)
        %Only one root
        z = myRoots(1);     
    else
        %Selection of Root which gives lower Gibb's energy
        tempVar = bb./bMix.*(myRoots(2)-myRoots(1)) - log((myRoots(2)-B_ofP) ...
             ./(myRoots(1)-B_ofP))- A_ofP./ (2.0*sqrt(2).*B_ofP).*(2.0*x*aij ... 
             ./aMix-bb./bMix).*log((myRoots(2)+2.414*B_ofP).*(myRoots(1)-0.414*B_ofP) ...
             ./((myRoots(2)-0.414*B_ofP).*(myRoots(1)+2.414*B_ofP)));
        dG = RR*TinK*(x*tempVar');
        if (dG > 0)
            z = myRoots(1);
        else
            z = myRoots(2);
        end
    end
    
%     v = z*RR*TinK/pres;    
%     denom = v*(v+bMix)+bMix*(v-bMix);
%     partialV = (RR*TinK/(v-bMix) + RR*TinK.*bb/((v-bMix)^2) ...
%      -2*x*aij/denom +(2*aMix*(v-bMix).*bb)/(denom*denom)) ... 
%      / (RR*TinK/((v-bMix)^2) -2*aMix*(v+bMix)/(denom*denom));
 
    logPhi = bb./bMix*(z-1.0) - log(z-B_ofP)-A_ofP./ ...
             (2.0*sqrt(2).*B_ofP).*(2.0*x*aij ... 
             ./aMix-bb./bMix).*log((z+2.414213562373095*B_ofP) ...
             ./(z-0.414213562373095*B_ofP));
         
         
    
%     if(bb>vMax || bb>v)
%         fprintf('Unphysical Result!!!\n');
%     end
end

