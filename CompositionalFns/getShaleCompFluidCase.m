
function [fluid, info] = getShaleCompFluidCase(name, varargin)
% Grab bag of different fluids for use in examples etc
    switch(lower(name))
        case 'lumped_1'
            % From SPE 79691 ex 5
            T_c = [189.515, 304.200, 387.607, 597.497, 698.515, 875.00];
            P_c = [45.2012, 72.90, 40.4196, 33.0150, 17.4525, 11.5372]*atm;
            mw = [16.1594, 44.01, 45.5725, 117.740, 248.827, 481.520]/1000;
            acc = [0.00854, 0.228, 0.16733, 0.38609, 0.80784, 1.23141];
            Z_c = [0.28981, 0.27055, 0.27588, 0.25668, 0.21967, 0.18250];
            
            V_c = Z_c.*8.314.*T_c./P_c;
            
            names = {'N2/CH4', 'CO2', 'C2-5', 'C6-13', 'C14-24', 'C25-80'};
            fluid = CompositionalFluid(names, T_c, P_c, V_c, acc, mw);
            
            bic = [0.11883, 0, 0, 0,0, 0;...
                   0.00070981, 0.15,0, 0, 0, 0; ...
                   0.00077754, 0.15,0, 0, 0, 0; ...
                   0.01, 0.15, 0, 0,0, 0; ...
                   0.011, 0.15, 0, 0,0, 0; ...
                   0.011, 0.15, 0, 0,0, 0];
               
            bic = bic + tril(bic, -1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0, 1, 0, 0, 0, 0], ...
                            'initial',   [0.463, 0.01640, 0.20520, 0.19108, 0.08113, 0.04319], ...
                            'pressure', 225*atm, ...
                            'temp', 387.45);
        case 'oil_1'
            %this fluid is used to test SlotDrill injector code
            names = {'Methane', 'CarbonDioxide', 'n-Decane'};
            
            info = makeInfo('injection', [0, 1, 0], ...
                            'initial',   [0.25, 0.25, 0.5], ...
                            'pressure', 5700*psia, ...
                            'temp', 352.594); %5700*psia
                        
            order = 1:3;
            info.injection = info.injection(order);
            info.initial = info.initial(order);
            fluid = TableCompositionalFluid(names(order));
            
            ChemPropsDB = ChemicalProps(names);
            bic = ChemPropsDB.bic;
            fluid = fluid.setBinaryInteraction(bic);
            case 'simple_1'
            %this fluid is used to test SlotDrill injector code
            names = {'Methane', 'n-Pentane', 'n-Decane'};
            
            info = makeInfo('injection', [1, 0, 0], ...
                            'initial',   [0.25, 0.25, 0.5], ...
                            'pressure', 5700*psia, ...
                            'temp', 352.594); %5700*psia
                        
            order = 1:3;
            info.injection = info.injection(order);
            info.initial = info.initial(order);
            fluid = TableCompositionalFluid(names(order));
            
            ChemPropsDB = ChemicalProps(names);
            bic = ChemPropsDB.bic;
            fluid = fluid.setBinaryInteraction(bic);
        case 'oil_1_modified'
            %this fluid is used to test SlotDrill injector code when we add
            %N2 to utilize GenericNatVars AD to calculate fugacity
            %derivatives
            names = {'Methane', 'CarbonDioxide', 'n-Decane','Nitrogen'};
            
            info = makeInfo('injection', [1,0, 0, 0], ...
                            'initial',   [0.25,0.25, 0.5, 0], ...
                            'pressure', 5700*psia, ...
                            'temp', 352.594); %5700*psia
                        
            order = 1:4;
            info.injection = info.injection(order);
            info.initial = info.initial(order);
            fluid = TableCompositionalFluid(names(order));
            
            ChemPropsDB = ChemicalProps(names);
            bic = ChemPropsDB.bic;
            fluid = fluid.setBinaryInteraction(bic);
            case 'oil_1_modified_2'
            %this fluid is used to test SlotDrill injector code when we add
            %N2 to utilize GenericNatVars AD to calculate fugacity
            %derivatives
            names = {'Nitrogen','Methane', 'CarbonDioxide', 'n-Decane'};
            
            info = makeInfo('injection', [1,0, 0, 0], ...
                            'initial',   [0,0.25, 0.25, 0.5], ...
                            'pressure', 5700*psia, ...
                            'temp', 352.594); %5700*psia
                        
            order = 1:4;
            info.injection = info.injection(order);
            info.initial = info.initial(order);
            fluid = TableCompositionalFluid(names(order));
            
            ChemPropsDB = ChemicalProps(names);
            bic = ChemPropsDB.bic;
            fluid = fluid.setBinaryInteraction(bic);
            case 'simple_1'
            %this fluid is used to test SlotDrill injector code
            names = {'Methane', 'n-Pentane', 'n-Decane'};
            
            info = makeInfo('injection', [1, 0, 0], ...
                            'initial',   [0.25, 0.25, 0.5], ...
                            'pressure', 5700*psia, ...
                            'temp', 352.594); %5700*psia
                        
            order = 1:3;
            info.injection = info.injection(order);
            info.initial = info.initial(order);
            fluid = TableCompositionalFluid(names(order));
            
            ChemPropsDB = ChemicalProps(names);
            bic = ChemPropsDB.bic;
            fluid = fluid.setBinaryInteraction(bic);
        case 'spe5'
            % Fifth SPE benchmark
            T_c = [343, 665.7, 913.4, 1111.8, 1270.0, 1380.0]*Rankine;
            P_c = [667.8, 616.3, 436.9, 304.0, 200.0, 162.0]*psia;
            mw = [16.040, 44.100, 86.180, 142.290, 206.000, 282.000]/1000;
            acc = [0.0130, 0.1524, 0.3007, 0.4885, 0.6500, 0.8500];
            Z_c = [0.290, 0.277, 0.264, 0.257, 0.245, 0.235];
            
            V_c = Z_c.*8.314.*T_c./P_c;
            
            names = {'C1', 'C3', 'C6', 'C10', 'C15', 'C20'};
            fluid = CompositionalFluid(names, T_c, P_c, V_c, acc, mw);
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(5, 1) = 0.05;
            bic(6, 1) = 0.05;
            bic(5, 2) = 0.005;
            bic(6, 2) = 0.005;
            
            bic = bic + tril(bic, -1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0.77, 0.20, 0.03, 0, 0, 0], ...
                            'initial',   [0.5, 0.03, 0.07, 0.20, 0.15, 0.05], ...
                            'pressure', 4000*psia, ...
                            'temp', 344.26);

        case 'fourmix'
            names = {'Methane', 'Nitrogen', 'n-Pentane', 'n-Decane'};
            
            info = makeInfo('injection', [0, 1, 0, 0], ...
                            'initial',   [0.1, 0.0, 0.45, 0.45], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 50);
                        
            fluid = TableCompositionalFluid(names);
            
%             ChemPropsDB = ChemicalProps(names);
%             bic = ChemPropsDB.bic;
%             fluid = fluid.setBinaryInteraction(bic);
            
        case '20components'
            names = {'Methane', 'Nitrogen', 'n-Pentane', 'n-Decane', 'CycloHexane', ...
                     'CarbonMonoxide', '1-Butene', 'Ammonia', 'HydrogenSulfide', 'IsoButane', ...
                     'n-Dodecane', 'n-Heptane', 'n-Hexane', 'n-Nonane', 'n-Octane', ...
                     'p-Xylene', 'n-Propane', 'n-Undecane', 'Oxygen', 'Water'};
            info = makeInfo('injection', (1:20)/20, ...
                            'initial',   ones(1, 20)/20, ...
                            'pressure', 300*barsa, ...
                            'temp', 273.15 + 50);
                        
            fluid = TableCompositionalFluid(names);

        case {'liquid_initial', 'vapor_initial'}
            names = {'Oxygen', 'Water'};
            if strcmpi(name, 'liquid_initial')
                inj =  [1, 0];
                init = [0, 1];
            else
                inj =  [0, 1];
                init = [1, 0];
            end            
            info = makeInfo('injection', inj, ...
                            'initial',   init, ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 30);
                        
            fluid = TableCompositionalFluid(names);
            fluid.names{2} = 'H2O';
            fluid.names{1} = 'O2';
        case 'watertracer'
            fluid = TableCompositionalFluid({'water', 'water'});
            info = makeInfo('injection', [1, 0], ...
                            'initial',   [0, 1], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 30);
            fluid.names = {'Tracer1', 'Tracer2'};
        case 'simple'
            names = {'Methane', 'n-Heptane', 'n-Decane'};
            
            info = makeInfo('injection', [0.1, 0.9, 0], ...
                            'initial',   [0.3, 0.1, 0.6], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 150);
                        
            order = 1:3;
            info.injection = info.injection(order);
            info.initial = info.initial(order);
            fluid = TableCompositionalFluid(names(order));
            
%             ChemPropsDB = ChemicalProps(names);
%             bic = ChemPropsDB.bic;
%             fluid = fluid.setBinaryInteraction(bic);
            
        case 'verysimple'
            names = {'CarbonDioxide', 'n-Decane'};
            fluid = TableCompositionalFluid(names);
            info = makeInfo('injection', [1, 0], ...
                            'initial',   [0.1, 0.9], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 150);
                        
%             ChemPropsDB = ChemicalProps(names);
%             bic = ChemPropsDB.bic;
%             fluid = fluid.setBinaryInteraction(bic);
        case 'onlydecane'
            names = {'n-Decane'};
            fluid = TableCompositionalFluid(names);
            info = makeInfo('injection', [1], ...
                            'initial',   [1], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 150);
                         
        case 'watermethane'
            names = {'Methane', 'Water'};
            fluid = TableCompositionalFluid(names);
            info = makeInfo('injection', [1, 0], ...
                            'initial',   [0.01, 0.99], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 150);  
                        
                        
        case 'onlymethane'
            names = {'Methane'};
            fluid = TableCompositionalFluid(names);
            info = makeInfo('injection', [1], ...
                            'initial',   [1], ...
                            'pressure', 5000*psia, ...
                            'temp', 366.48);
        case 'c1c2'
            names = {'Methane','Ethane'};
            fluid = TableCompositionalFluid(names);
            info = makeInfo('injection', [0.999 0.001], ...
                            'initial',   [0.999 0.001], ...
                            'pressure', 5000*psia, ...
                            'temp', 366.48); 
                     
        case 'barnett3comps'
            %OMO: Edit added a couple of hydrocarbon mixtures
            % 3 component light-gas mixture
            names = {'Methane', 'Ethane', 'n-Propane'};
            %Pressure in Pa in first column, rho_sL in kg/m^3 in 2nd column
            isotherm = [10769611, 2.99;  5591648, 4.86; 5819175, 9.56];
            fluid = SorbedCompositionalFluid(names,isotherm');   
            
            Di=[2.8,2.5,1.9]*10^-7; %SPE-175074-MS
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 2) = 0.005;
            bic(1, 3) = 0.01;
            bic(2, 3) = 0.005;
            
            bic = bic + triu(bic, 1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0, 0, 1.0], ...
                            'initial',   [0.991,0.0088,0.0002], ...
                            'pressure', 5000*psia, ...
                            'temp', 355.84); 
                        
        case 'barnett3comps_modified_2'
            %Hassan: Edit added CO2 to the last to work with GenericNatVars
            % 3 component light-gas mixture + imaginary CO2
            names = {'Methane', 'Ethane', 'n-Propane','CarbonDioxide'};
            %Pressure in Pa in first column, rho_sL in kg/m^3 in 2nd column
            isotherm = [ 10769611, 2.99;  5591648, 4.86; 5819175, 9.56;10769611, 2.99];
            fluid = SorbedCompositionalFluid(names,isotherm');  
           
           
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 2) = 0.005;
            bic(1, 3) = 0.01;
            bic(1, 4) = 0;
            bic(2, 3) = 0.01;
            bic(2, 4) = 0;
            bic(3, 4) = 0;
           
            bic = bic + triu(bic, 1)';
           
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [1, 0,0 , 0], ...
                            'initial',   [0.991,0.0088,0.0002,0], ...
                            'pressure', 5000*psia, ...
                            'temp', 355.84);
                        
        case 'barnett3comps_nat_vars'
            %Hassan: Edit added CO2 to the last to work with GenericNatVars
            % 3 component light-gas mixture + imaginary CO2
            names = {'Methane', 'Ethane', 'n-Propane','CarbonDioxide'};
            %Pressure in Pa in first column, rho_sL in kg/m^3 in 2nd column
            isotherm = [ 10769611, 2.99;  5591648, 4.86; 5819175, 9.56;10769611, 2.99];
            fluid = SorbedCompositionalFluid(names,isotherm');  
           
           
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 2) = 0.005;
            bic(1, 3) = 0.01;
            bic(1, 4) = 0;
            bic(2, 3) = 0.01;
            bic(2, 4) = 0;
            bic(3, 4) = 0;
           
            bic = bic + triu(bic, 1)';
           
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [1, 0,0 , 0], ...
                            'initial',   [0.991,0.0088,0.0002,0], ...
                            'pressure', 5000*psia, ...
                            'temp', 355.84);
         case 'barnett3comps_modified'
            %Hassan: Edit added CO2 to calculate fugacity derivatives
            % 3 component light-gas mixture + imaginary CO2
            names = {'CarbonDioxide','Methane', 'Ethane', 'n-Propane'};
            %Pressure in Pa in first column, rho_sL in kg/m^3 in 2nd column
            isotherm = [10769611, 2.99; 10769611, 2.99;  5591648, 4.86; 5819175, 9.56];
            fluid = SorbedCompositionalFluid(names,isotherm');   
            
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 2) = 0;
            bic(1, 3) = 0;
            bic(1, 4) = 0;
            bic(2, 3) = 0.005;
            bic(2, 4) = 0.01;
            bic(3, 4) = 0.005;
            
            bic = bic + triu(bic, 1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0, 0, 0, 1.0], ...
                            'initial',   [0, 0.991,0.0088,0.0002], ...
                            'pressure', 5000*psia, ...
                            'temp', 355.84); 
        case 'c1c2c3'
            %OMO: Edit added a couple of hydrocarbon mixtures
            % 3 component light-gas mixture
            names = {'Methane', 'Ethane', 'n-Propane'};
            fluid = TableCompositionalFluid(names);
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 2) = 0.005;
            bic(1, 3) = 0.01;
            bic(2, 3) = 0.005;
            
            bic = bic + triu(bic, 1)';
%             ChemPropsDB = ChemicalProps(names);
%             bic = ChemPropsDB.bic;
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0, 0, 1.0], ...
                            'initial',   [0.991,0.0088,0.0002], ...
                            'pressure', 5466.482220000024*psia, ...
                            'temp', 355.84); 
                          
        case 'bakken'
            % Bakken shale composition from Nojabaei dissertation 2015
            T_c = [335.336	549.969	665.97	759.208	875.479	1053.25	1332.095	1844.491].*Rankine();
            P_c = [655.02	721.99	615.76	546.46	461.29	363.34	249.61	190.12]*psia;
            mw = [16.535	30.433	44.097	58.124	78.295	120.562	220.716	443.518]/1000;
            acc = [0.0102	0.1028	0.152	0.1894	0.2684	0.4291	0.7203	1.0159];
            
            V_c = [1.58	2.34	3.25	4.11	5.39	8.81	15.19	36].*0.06242795996802./1000;
            Z_c = P_c.* V_c./8.314./T_c;
            Parachor = [74.8	107.7	151.9	189.6	250.2	350.2	590	1216.8];
            
            names = {'C1','C2','C3','C4','C5-C6','C7-C12','C13-C21','C22-C80'};
            fluid = CompositionalFluid(names, T_c, P_c, V_c, acc, mw);
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 1:8) = [0	0.005	0.0035	0.0035	0.0037	0.0033	0.0033	0.0033];
            bic(2, 1:8) = [0	0	0.0031	0.0031	0.0031	0.0026	0.0026	0.0026];

            bic = bic + triu(bic, 1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [1	0 0 0 0 0 0 0 ], ...
                            'initial',   [0.36736	0.14885	0.09334	0.05751	0.06406	0.15854	0.0733	0.03704], ...
                            'pressure', 5700*psia, ...
                            'temp', 352.594);
        case 'eagleford'
            % HR-Edit-EagleFord shale composition from 
            T_c = [547.56 227.16 343.08 494.532 789.624	1332.522].*Rankine();
            P_c = [1069.865	492.3143 667.1961 536.4021 368.5744 257.9139]*psia;
            mw = [44.01	28.01 16.04	52.02 103.01 267.15]/1000;
            acc = [0.225 0.04 0.008	0.1723	0.2839	0.6716];

            V_c = [9.4e-5 8.95e-5 9.9e-5 0.0002293 0.0003943 0.000887]; %m3/mol
            Z_c = P_c.* V_c./8.314./T_c;
            Parachor = [78	41	77	171.07	297.42	661.45];
            
            names = {'CO2','N2','C1','C2-C5','C6-C10','C11+'};
            fluid = CompositionalFluid(names, T_c, P_c, V_c, acc, mw);
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 1:6) = [0    0.02    0.103	0.1299  0.15    0.15];
            bic(2, 1:6) = [0	0	0.031	0.082	0.12	0.12];
            bic(3, 1:6) = [0	0	0	0.0174	0.0462	0.111];
            bic(4, 1:6) = [0	0	0	0	0.0073	0.0444];
            bic(5, 1:6) = [0	0	0	0	0	0.0162];


            bic = bic + triu(bic, 1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0	0 1 0 0 0], ...
                            'initial',   [0.01183	0.00161	0.11541	0.26438	0.38089	0.22588], ...
                            'pressure', 8125*psia, ...
                            'temp', 405.372); % The Effect of Temperature on Flowback Data Analysis in Shale Gas Reservoirs(yang et al); 

       case 'bakken_light'
            % Bakken shale composition from Yu et al. 2015
            % This is a simplified Bakken fluid and it is used in SlotDril
            % EOR Injection Code
            
            T_c = [304.2 126.2 190.6 363.3 511.56 579.34 788.74];
            P_c = [72.8 33.5 45.4 42.54 33.76 30.91 21.58]*atm;
            mw = [44.01 28.01 16.04 42.82 83.74 105.91 200]./1000;
            acc = [0.2250 0.0400 0.008 0.1432 0.2474 0.2861 0.6869];
            
            V_c = [0.0940 0.0895 0.0990 0.1970 0.3338 0.4062 0.9208]./1000;
            Z_c = P_c.* V_c./8.314./T_c;
            Parachor = [78 41 77 145.2 250 306 686.3];
            
            names = {'CO2','N2','C1','C2-C4','C5-C7','C8-C9','C10+'};
            fluid = CompositionalFluid(names, T_c, P_c, V_c, acc, mw);
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 1:7) = [0 -0.02 0.103 0.1327 0.1413 0.15 0.15];
            bic(2, 1:7) = [-0.02 0 0.013 0.0784 0.1113 0.12 0.12];
            bic(3, 1:7) = [0.103 0.0130 0	0.0078	0.0242	0.0324	0.0779];
            bic(4, 1:7) = [0.1327 0.0784 0.0078	0	0.0046	0.0087	0.0384];
            bic(5, 1:7) = [0.1413 0.1113 0.0242	0.0046	0	0.0006	0.0169];
            bic(6, 1:7) = [0.15 0.12 0.0324	0.0087	0.0006	0	0.0111];
            bic(7, 1:7) = [0.15 0.12 0.0779	0.0384	0.0169	0.0111	0];
            
            fluid = fluid.setBinaryInteraction(bic);
%             'initial',   [0.2506 0.22 0.2 0.13 0.1994], ...
            info = makeInfo('injection', [0 0 1 0 0 0 0], ...
                            'initial',   [0.0002, 0.0004, 0.25, 0.22, 0.2, 0.13, 0.1994], ...
                            'pressure', 4000*psia, ...
                            'temp', 388.706);       %1500*psia          
       case 'bakken_light_dummy'
            % Bakken shale composition from Yu et al. 2015
            % This is a simplified Bakken fluid and it is used to utilize
            % GenericNatVars approach
            
            T_c = [304.2 126.2 190.6 363.3 511.56 579.34 788.74 100];
            P_c = [72.8 33.5 45.4 42.54 33.76 30.91 21.58 10]*atm;
            mw = [44.01 28.01 16.04 42.82 83.74 105.91 200 50]./1000;
            acc = [0.2250 0.0400 0.008 0.1432 0.2474 0.2861 0.6869 0.65];
            
            V_c = [0.0940 0.0895 0.0990 0.1970 0.3338 0.4062 0.9208 0.6]./1000;
           Z_c = P_c.* V_c./8.314./T_c;
            Parachor = [78 41 77 145.2 250 306 686.3 80];
            
            names = {'CO2','N2','C1','C2-C4','C5-C7','C8-C9','C10+','dummy'};
            fluid = CompositionalFluid(names, T_c, P_c, V_c, acc, mw);
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 1:8) = [0 -0.02 0.1030 0.1327 0.1413 0.15 0.15 0.1];
            bic(2, 1:8) = [-0.02 0 0.013 0.0784 0.1113 0.12 0.12 0.2];
            bic(3, 1:8) = [0.103 0.0130 0	0.0078	0.0242	0.0324	0.0779 0.3];
            bic(4, 1:8) = [0.1327 0.0784 0.0078	0	0.0046	0.0087	0.0384 0.4];
            bic(5, 1:8) = [0.1413 0.1113 0.0242	0.0046	0	0.0006	0.0169 0.5];
            bic(6, 1:8) = [0.15 0.12 0.0324	0.0087	0.0006	0	0.0111 0.6];
            bic(7, 1:8) = [0.15 0.12 0.0779	0.0384	0.0169	0.0111	0 0.7];
            bic(8, 1:8) = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0];
            
            fluid = fluid.setBinaryInteraction(bic);
%             'initial',   [0.2506 0.22 0.2 0.13 0.1994], ...
            info = makeInfo('injection', [1 0 0 0 0 0 0 0], ...
                            'initial',   [0.0002, 0.0004, 0.25, 0.22, 0.2, 0.13, 0.1994 0], ...
                            'pressure', 8000*psia, ...
                            'temp', 388.706);       %1500*psia         
                        
        otherwise
            error('Unknown case')
    end
end

function info = makeInfo(varargin)
    info = struct('injection', [], 'initial', [],...
                  'pressure',  [], 'temp',    []);
    info = merge_options(info, varargin{:});
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.
This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}