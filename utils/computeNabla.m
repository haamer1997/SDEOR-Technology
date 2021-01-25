function dd = computeNabla(G, rock) %denominator
    opt = struct('verbose'        , false, ...
                'grdecl'         , []   , ...
                'K_system'       , 'xyz', ...
                'cellCenters'    , []   , ...
                'fixNegative'    , true, ...
                'cellFaceCenters', []   );
    [C, N, cellNo] = cell_geometry(G, opt);
    [K, i, j] = permTensor(rock, G.griddim);

    assert (size(K,1) == G.cells.num, ...
          ['Permeability must be defined in active cells only.\n', ...
           'Got %d tensors, expected %d (== number of cells).'],   ...
           size(K,1), G.cells.num);

    % Compute T = C'*K*N / C'*C. Loop-based to limit memory use.
    unitNorm = N ./ abs(N); 
    unitNorm(isnan(unitNorm)) = 0; % get
    N = unitNorm;
    unitK = K ./ K; 
    unitK(isnan(unitK)) = 0; % get
    dd = zeros(size(cellNo));
    for k = 1 : size(i, 2)
      dd = dd + (C(:, i(k)) .* unitK(cellNo, k) .* N(:, j(k)));
    end

    dd = dd ./ sum(C .* C, 2);
    dd  = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ dd, [G.faces.num, 1]); 
    
    if isfield(G,'nnc')
        tol=1e-5;
        FracGrid=G.FracGrid;
        for i=1:numel(fieldnames(FracGrid)) % loop through every fracture
            fieldname=['Frac',num2str(i)];
            F=FracGrid.(fieldname);

            % Detect NNCs
            [cellsi,CIi,~,typei]=fracmatnnc(F,G.Matrix,tol);

            % Compute NNC Transmissibilities
%             w1 = pv(cellsi(:,1))./G_global.rock.perm(cellsi(:,1));
%             w2 = pv(cellsi(:,2))./G_global.rock.perm(cellsi(:,2));
%             wt = pv(cellsi(:,1))+pv(cellsi(:,2));
            % No weighting by PV:
            w1 = 1; %1./G.rock.perm(G.nnc.cells(:,1)); %Assume 1
            w2 = 1; %1./G.rock.perm(G.nnc.cells(:,2));
            wt = 1;
            Ti = CIi.*(wt./(w1+w2));

            % Adjust Ti for surface area. Ti was previously calculated using
            % frac-mat intersection area. For fractures within a cell, surface area
            % is 2*intersection. For fractures on a cell face, surface area is
            % equal to intersection area. Refer to Tene et al. (2017) pEDFM.
%             Ti = adjustT(Ti, cellsi); %do I need it to calculate nabla?

            % Append data to existing nnc information
            dd=[dd;Ti];  
        end
    end
end

function [C, N, cellNo] = cell_geometry(G, opt)
   % Vectors from cell centroids to face centroids
   cellNo = gridCellNo(G);

   if ~isempty(opt.cellCenters)
      C = opt.cellCenters;
   else
      C = G.cells.centroids;
   end

   if ~isempty(opt.cellFaceCenters)
      C = opt.cellFaceCenters - C(cellNo,:);
   else
      C = G.faces.centroids(G.cells.faces(:,1), :) - C(cellNo,:);
   end

   % Outward-pointing normal vectors
   cf  = G.cells.faces(:,1);
   sgn = 2*(cellNo == G.faces.neighbors(cf, 1)) - 1;
   N   = bsxfun(@times, sgn, G.faces.normals(cf, :));
end