% Class for the Lancet Countdown Europe Mobility Indicator 2022 
% Henrik Sjödin (2022), henrik.sjodin@umu.se

% ToDo: 
% (1*) In builder; create export- and import location already in builder, and store in dataTable that is passed to object MobilityIndicator.
% (2*) In MobilityIndicator; create function that computes per-location exports and imports for year y. Store as properties in object MobilityIndicator. 
% (3*) In MobilityIndicator; rewrite worldMap() and drawPolygon() accordingly;
% (4) Create special purpose function; write functions that plot time-trends. 
% (5) Sort essential Matlab- and data files, tidy in a CORE-folder.
% (6*) Exclude European countires from coutry data as nuts data replace
% this.
% (7) Commit v1.0 to git

classdef MIBuilder_LCDEurope
    properties (Hidden)
        year
        travellerProportion (1,:)
        incidenceThreshold
        countryCodesConversionTable
        worldCountryShapes
        nutsShapes
        nutsData

    end
    methods
        function this = MIBuilder_LCDEurope(year,travellerProportion, incidenceThreshold)
            this.year = year;
            this.travellerProportion = travellerProportion;
            this.incidenceThreshold = incidenceThreshold;
            this.countryCodesConversionTable = readtable('Data/country_codes.txt');
            this.worldCountryShapes = shaperead('CNTR_RG_03M_2014.shp', 'UseGeoCoords', true);
            
            nutsData = readtable('Data/NUTS_region_R0_LTS_population_1950_2021_2022-02-26.csv');
            nutsData = nutsData(nutsData.year == this.year,:);
            nid = arrayfun(@(x) string(x), nutsData.NUTS_ID, 'UniformOutput', false);
            cid = arrayfun(@(x) string(x), nutsData.CNTR_CODE, 'UniformOutput', false);
            this.nutsData = nutsData([nid{:}]'~=[cid{:}]',:); % Exclude countries (if not equal to)
            clear nid cid nutsData
            
            nutsShapes = shaperead('Euro_Nuts_regions.shp', 'UseGeoCoords', true);
            nid = arrayfun(@(x) string(x.NUTS_ID), nutsShapes, 'UniformOutput', false);
            cid = arrayfun(@(x) string(x.CNTR_CODE), nutsShapes, 'UniformOutput', false);
            this.nutsShapes = nutsShapes([nid{:}]'~=[cid{:}]'); % Exclude countries (if not equal to)
            
            
            clear nid cid nutsShapes
        end
        
        function buildObj = build(this)
            dataTable = this.prepareDataTable(readtable('Data/Country_Population_GBD_dengue_incidence_rate_1990_2020_2022-02-21_edited.xlsx'));
%             destinationLocationTable = readtable('european_countires.txt');
            for i = 1:size(dataTable,1)
                travProp(1,i) = this.travellerProportion.getTravellerProportion(dataTable(i,:).Country,this.year);
            end
            infectedTravellerProportion = travProp .* dataTable.infectedProportion';
            buildObj = MobilityIndicator(dataTable, this.year, infectedTravellerProportion, this.incidenceThreshold, this.worldCountryShapes, this.nutsShapes);
            clear dataTable destinationLocationTable worldCountryShapes
        end
    end
    
    methods
        
        function dataTable = prepareDataTable(this, dataTable) % This function may need to be slightly reformulated depending on input data formats!!!
%             dataTable = renamevars(dataTable,["x","y"],["lon","lat"]); % specifically for first tentative datatable?
            dataTable = renamevars(dataTable,["ISO3"],["CODE"]);
            dataTable = dataTable(dataTable.year==this.year,:);
            ecd = readtable('Data/european_countires.txt');
            nonEuropeanCountriesL = ~ismember(dataTable.Country,ecd.name);
            dataTable = dataTable(nonEuropeanCountriesL,:); % Exclude european countries (as these locations are covered by nuts-regions instead)
            cid=(1:size(dataTable,1))';
            CID = table(cid);
            dataTable = [CID dataTable]; clear T0
            zCL = dataTable.Population==0; % specifically for first tentative datatable?
            if isempty(dataTable)
                error('The input year is not represented in the input data table. use a different year or a different data.')
            elseif sum(zCL)>0
                dataTable = dataTable(~zCL,:); warning('%i countries with zero population size exist, and %i are removed',sum(zCL),sum(zCL))
            end
%             this.dataTable = dataTable;

            
            cShapes = this.worldCountryShapes;
            for i = 1:size(dataTable,1)
                iso2 = this.iso3ToIso2( dataTable(i,:).CODE );%this.iso3ToIso2( dataTable(i,:).ISO3 );
                for j = 1:numel(cShapes)
                    if strcmp(cShapes(j).CNTR_ID,iso2{1})
%                         disp('prepareDataTable(): Found!')
                        COUNTRYSHP(i,1) = cShapes(j); 
                        NUTSSHP(i,1) = structfun(@(x) [], this.nutsShapes(1), 'UniformOutput', false);
                        cShapes(j) = [];
                        break
                    elseif j == numel(cShapes)
                        %warning(['No shapefile-match for %s.'],
                        %dataTable(i,:).Country{1} ) %%%%%% Write this to log-file instead!
                        COUNTRYSHP(i,1) = structfun(@(x) [], COUNTRYSHP(i-1,1), 'UniformOutput', false);
                        NUTSSHP(i,1) = structfun(@(x) [], this.nutsShapes(1), 'UniformOutput', false);
                        
                    end
                end
            end
            NUTS = zeros(size(dataTable,1),1);
            R0 = zeros(size(dataTable,1),1);
            LTS = R0;
            
            dataTable = [dataTable table(COUNTRYSHP) table(NUTSSHP) table(NUTS) table(R0) table(LTS) ];
            dataTable = this.addNUTS3ToDataTable(dataTable);
            logicalExportLocations = dataTable.dengue_incidence_rate > this.incidenceThreshold; % set export location
            logicalImportLocations = ( (dataTable.NUTS == 1) & (dataTable.LTS >= 1.0)); % set import location
            infectedProportion = dataTable.dengue_incidence_rate * 10 / (365 * 100000); % incidence is in the unit number of cases per 100k and per year. MAKE PER DAY, AND ASSUME A 10 DAY WINDOW. Total population-size cancels
            dataTable = [dataTable table(logicalExportLocations) table(logicalImportLocations) table(infectedProportion)];
        end
        
        function dataTable = addNUTS3ToDataTable(this, dataTable)
            cid = (1:size(this.nutsData,1))' + dataTable.cid(end);
            CODE = this.nutsData.NUTS_ID;
            year = this.nutsData.year;
            Country = this.nutsData.NUTS_NAME;
            lon = this.nutsData.x;
            lat = this.nutsData.y;
            Population = this.nutsData.population;
            dengue_incidence_rate = zeros(size(this.nutsData,1),1);
            
            
            
            T = table(cid,CODE,year,Country,lon,lat,Population,dengue_incidence_rate);
            
            nShapes = this.nutsShapes;
            
            for i = 1:size(T,1)
%                 iso2 = this.iso3ToIso2( dataTable(i,:).CODE );;%this.iso3ToIso2( dataTable(i,:).ISO3 );
                for j = 1:numel(nShapes)

                    if strcmp(nShapes(j).NUTS_ID,CODE(i))
                        COUNTRYSHP(i,1) = structfun(@(x) [], this.worldCountryShapes(1), 'UniformOutput', false);
                        NUTSSHP(i,1) = nShapes(j); 
                        nShapes(j) = [];
                        break
                    elseif j == numel(nShapes)
                        %warning(['No shapefile-match for %s.'], T(i,:).Country{1} ) %%%%%% Write this to log-file instead!
                        COUNTRYSHP(i,1) = structfun(@(x) [], this.worldCountryShapes(1), 'UniformOutput', false);
                        NUTSSHP(i,1) = structfun(@(x) [], this.nutsShapes(1), 'UniformOutput', false);
                        
                    end
                end
            end
            
            NUTS = ones(size(T,1),1);
            R0 = this.nutsData.R0;
            LTS = this.nutsData.LTS;
            T = [T table(COUNTRYSHP) table(NUTSSHP) table(NUTS) table(R0) table(LTS)];

            dataTable = [dataTable; T];
            
        end
        
        function iso = iso3ToIso2(this, iso3)
            iso = this.countryCodesConversionTable(strcmp( this.countryCodesConversionTable.alpha_3, iso3),:).alpha_2;
        end
    end
end

