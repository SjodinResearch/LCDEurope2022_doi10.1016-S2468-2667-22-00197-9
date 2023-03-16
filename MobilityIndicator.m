% Class for the Lancet Countdown Mobility Indicator 2022 
% Henrik Sjödin (2022), henrik.sjodin@umu.se

classdef MobilityIndicator
    properties
        dataTable table
        year
        travellerProportion (1,:) double
        
        incidenceThreshold
        worldCountryShapes
        nuts3Shapes
        sMatrix
        length_of_visitation
    end
    methods 
        function this = MobilityIndicator(dataTable, year,travellerProportion, incidenceThreshold, worldCountryShapes, nuts3Shapes, length_of_visitation)
            arguments
                dataTable table
                year double
                travellerProportion double {mustBeGreaterThanOrEqual(travellerProportion,0),mustBeLessThanOrEqual(travellerProportion,1)}
                incidenceThreshold double {mustBeGreaterThanOrEqual(incidenceThreshold,0)}
%                 destinationLocationTable table
                worldCountryShapes struct
                nuts3Shapes struct
                length_of_visitation double
            end
            this.year = year; clear year
            this.travellerProportion = travellerProportion; clear travellerProportion
            this.incidenceThreshold = incidenceThreshold; clear incidenceThreshold
            this.dataTable = dataTable;
            
            this.worldCountryShapes = worldCountryShapes;
            this.nuts3Shapes = nuts3Shapes;
            this.length_of_visitation = length_of_visitation;
            clear dataTable
        end
    end
    
    methods

        function d = distanceMatrix(this)
            % Returns "geodesic" distance in km from lat long ccordinates
            % Distance using matrices - efficient calculation
            ln_mat = [this.dataTable.lat, this.dataTable.lon];
            n = length(ln_mat);
            DM = [repelem(ln_mat, n, 1) repmat(ln_mat, n, 1)];
            dm = deg2km(distance(DM(:, 1), DM(:, 2), DM(:, 3), DM(:, 4)));
            d = reshape(dm, [n n]); % reshape
            clear ln_mat n DM dm
        end
%         function s = sMatrixSlow(this) % Obsolete
%             tic;
%             d = this.distanceMatrix;
%             n = this.dataTable.Population;%populationSizes;
%             s = zeros(size(d));
%             for i = 1:size(d, 1)
% 
%                 for j = 1:size(d, 1)
% 
%                     if i == j
%                         s(i, j) = 0;
%                     else
%                         indx = find(d(i, :) < d(i, j) & d(i, :) ~= d(i, i));
% 
%                         if (isempty(indx) == 1)
%                             s(i, j) = 0;
%                             %T(i,j)= (Ni*Nj)/((Ni+s(i,j))*(Ni+Nj+s(i,j)));
%                         else
%                             s(i, j) = sum(n(indx'));
%                             %T(i,j)= (Ni*Nj)/((Ni+s(i,j))*(Ni+Nj+s(i,j)));
%                         end
% 
%                     end
% 
%                 end
% 
%             end
%             toc
% 
%         end
        
        function s = getsMatrix(this)
            if isempty(this.sMatrix)
                d = this.distanceMatrix;
                n = this.dataTable.Population;%populationSizes;
                s = zeros(size(d));
                parfor i = 1:size(d,1)
                    [xi,xj] = meshgrid(d(i,:));
                    s(i,:) = sum(n.*(xi>xj),1);
                end
            else
                s = this.sMatrix;
            end
        end
        function this = store_sMatrix(this)
            this.sMatrix = this.getsMatrix;
        end
        function [N, ni, nj, s] = prepareVectirizedVariablesForRAdiationModel(this)
            % N is population matrix
            % X is total number of individuals that potentially could travel
            N = this.dataTable.Population;%this.populationSizes;
            if isempty(this.sMatrix)
                error(' ** sMatrix must first be computed and made available by running obj = obj.store_sMatrix.')
%                 s = this.sMatrix;
            else
                s = this.sMatrix;
            end
            [ni, nj] = ndgrid(N, N);
        end
        function p = travelProbMat(this)%(N, s)%radiation_model_reallyvectorized(X,N,s)
            [N, ni, nj, s] = this.prepareVectirizedVariablesForRAdiationModel;
            p = (ni .* nj) ./ ((ni + s) .* (ni + nj + s));
            p(1:numel(N) + 1:end) = 0; % Diagonal?
        end
        function p = travelRateMat(this,travellerProportion)%(N, s)%radiation_model_reallyvectorized(X,N,s)
            [N, ni, nj, s] = this.prepareVectirizedVariablesForRAdiationModel;
            p = travellerProportion * ni .* (ni .* nj) ./ ((ni + s) .* (ni + nj + s));%this.travellerProportion * ni .* (ni .* nj) ./ ((ni + s) .* (ni + nj + s));
            p(1:numel(N) + 1:end) = 0; % Diagonal?
        end
        
    end
    
    methods

        function r = travelRatesIntoDestinationsFromSourceCountires(this, travellerProportion)
            r0 = this.travelRateMat(travellerProportion); %this.travelRateMat
%             r = r(this.destLogi,this.sourceLogi);
            r = r0(this.dataTable.logicalImportLocations,this.dataTable.logicalExportLocations);
            clear r0
        end
        function [exports,imports] = perLocationFlowRate(this, travellerProportion)
            
            r = this.travelRatesIntoDestinationsFromSourceCountires(travellerProportion);
            exports = sum(r,1)';
            imports = sum(r,2);
            
        end

        function r = travelRatesofVisitorsIntoDestinationsFromSourceCountires(this, travellerProportion)
            r0 = this.travelRateMat(travellerProportion); %this.travelRateMat
%             r = r(this.destLogi,this.sourceLogi);
            r = r0(this.dataTable.logicalExportLocations, this.dataTable.logicalImportLocations);
            clear r0
        end
        function imports = importsFromVistors(this)
            rv = = travelRatesofVisitorsIntoDestinationsFromSourceCountires(this, travellerProportion);
            imports = rv .* this.length_of_visitation .* "daily infetion rate in visited country"
        end
%         function v = importsPerMonthTimesLTS(this)
%             [~,imports] = perLocationFlowRate(this);
%             v = ( imports / 12 ) .* this.dataTable(this.dataTable.logicalImportLocations).LTS;
%         end
    end
    
    %%% MAPPING -----------------------------------------------------------
    
    methods
        function worldMap(this,fnum)
            ax = this.drawMap(fnum, 'World');
            setm(ax,'MapLatLimit',[-60 75])
            geoshow(ax,'landareas.shp', 'FaceColor',[0.9 0.98 0.9])
            
            if any(this.dataTable.logicalExportLocations) & any(this.dataTable.logicalImportLocations)
                [exp,imp] = this.perLocationFlowRate(this.travellerProportion);

                cmape = flipud(copper(128));%parula(numel(exp));
                cmapi = flipud(copper(128));

                eSHP = this.dataTable(this.dataTable.logicalExportLocations,:).COUNTRYSHP;
                arrayfun(@(SHP,t) this.drawPolygon(ax, SHP, this.scalar2Color(t,[min(exp) max(exp)],cmape)) , eSHP, exp);

                iSHP = this.dataTable(this.dataTable.logicalImportLocations,:).NUTSSHP;
                arrayfun(@(SHP,t) this.drawPolygon(ax, SHP, this.scalar2Color(t,[min(imp) max(imp)],cmapi)) , iSHP, imp)  ;
            else
                disp('No risk related travels.')
            end
            text(0.15,0.2,0,['Year: ' num2str(this.year)],'Units','normalized')
            colormap(cmape)
            colorbar(ax,'Ticks',[0,1],...
         'TickLabels',{'Low','High'})
        end
        
        function europeMap(this,fnum)
            ax = this.drawMap(fnum, 'Europe');
%             setm(ax,'MapProjection','eqdcylin','MapLatLimit',[40 70])
            geoshow(ax,'landareas.shp', 'FaceColor',[0.9 0.98 0.9])
            if any(this.dataTable.logicalExportLocations) & any(this.dataTable.logicalImportLocations)
                [~,imp] = this.perLocationFlowRate(this.travellerProportion);
                impLTS = (imp / 12) .* this.dataTable(this.dataTable.logicalImportLocations,:).LTS;
%                 cmape = flipud(copper(64));%parula(numel(exp));
                cmapi = flipud(copper(128));

%                 eSHP = this.dataTable(this.dataTable.logicalExportLocations,:).COUNTRYSHP;
%                 arrayfun(@(SHP,t) this.drawPolygon(ax, SHP, this.scalar2Color(t,[min(exp) max(exp)],cmape)) , eSHP, exp);

                iSHP = this.dataTable(this.dataTable.logicalImportLocations,:).NUTSSHP;
                arrayfun(@(SHP,t) this.drawPolygon(ax, SHP, this.scalar2Color(t,[min(impLTS) max(impLTS)],cmapi)) , iSHP, impLTS)  ;
            else
                disp('No risk related travels.')
            end
            colormap(cmapi)
            colorbar(ax,'Ticks',[0,1],...
         'TickLabels',{'Low','High'})
            text(0.1,0.15,0,['Year: ' num2str(this.year)],'Units','normalized')
        end

        function europeMapRiskRegions(this, fnum)
            cmap = flipud(copper(128));
            ax = this.drawMap(fnum, 'Europe');
            geoshow(ax,'landareas.shp', 'FaceColor',[0.9 0.98 0.9])
            iSHP = this.dataTable(this.dataTable.logicalImportLocations,:).NUTSSHP;
            LTS = this.dataTable(this.dataTable.logicalImportLocations,:).LTS;
            arrayfun(@(SHP,t) this.drawPolygon(ax, SHP, this.scalar2Color(t,[min(LTS) max(LTS)],cmap)) , iSHP, LTS);
           
            colormap(cmap)
            c = colorbar(ax,'Ticks',[0,1],...
         'TickLabels',{num2str(min(LTS)),num2str(max(LTS))})
            c.Label.String = 'Length of transmission season (months)';
            text(0.1,0.15,0,['Year: ' num2str(this.year)],'Units','normalized')
        end
        
        function ax = drawMap(this, fnum, region)
            figure(fnum); %clf(fnum)
            disp(['Drawing ' region ' map for year ' num2str(this.year) ', in figure ' num2str(fnum) ' ...'])
            ax = worldmap( region );
            setm(ax,'MapProjection','eqdcylin')
            mlabel off; plabel off;
        end
        
        function drawPolygon(this, ax, SHP, color)
%             disp(['Drawing polygon.'])
            geoshow(ax,SHP.Lat,SHP.Lon,'DisplayType','Polygon', 'FaceColor',color);
        end
        
        function c = scalar2Color(this, scalar, extremValues, RGBmatrix)
            ev =[floor(extremValues(1)) ceil(extremValues(2))];
            x = (ev(1):(range(ev)/(size(RGBmatrix,1) - 1)):ev(2))';
            c = interp1(x,RGBmatrix,scalar);
        end
        
    end
end