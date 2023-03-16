classdef TravellerProportion
    properties
        tProp
        countriesNames string
        intercept
        slope
    end
    methods
        function this = TravellerProportion 

            UNWTO = readtable('Data/UNWTO/UNWTO_Travel_1995_2011.xlsx');

            x = unique(UNWTO.dyear);
            % c = unique(UNWTO.name_origin);
            
            
            
            
%             m=MIBuilder_LCDEurope(1990,0.01, 500).build;
%             mT = m.dataTable;%(m.dataTable.NUTS == 0,:);
%             c = unique(mT(mT.NUTS==0,:).Country);
            T = readtable('Data/Country_Population_GBD_dengue_incidence_rate_1990_2020_2022-02-21_edited.xlsx');
               c = unique(T.Country);
            arrivals = zeros(numel(x),numel(c));
            pop = arrivals;
            % tProp = arrivals;
            match = logical(arrivals);
            
            k = 0;
            for i = 1:numel(x)
                t=UNWTO(UNWTO.dyear==x(i),:);
%                 m=MIBuilder_LCDEurope(x(i),0.01, 500).build;
                
                mT = T(T.year == x(i),:);%m.dataTable;%(m.dataTable.NUTS == 0,:);
            %     c = unique(mT(mT.NUTS==0,:).Country);
                for j = 1:numel(c)
                    arrivals(i,j) = sum(t(strcmp(t.name_origin,c(j)),:).arrivals);
                    I = find(strcmp(mT.Country,c(j)));
                    if ~isempty(I)
                        if numel(I) >= 2 
                            I = I(1);
                        end
                        pop(i,j) = mT(I,:).Population;
                        if i == 1
                            match(j) = logical(1);
                            nm = c(j);
                            this.countriesNames(j) = string(nm{1});
                        end
                    end
                end
            end
            
            for i = 1:numel(x)
                k = 0;
                for j = 1:numel(c)
                    if match(j)
                        k=k+1;
                        this.tProp(i,k) = arrivals(i,j) / pop(i,j);
                        

                    end
                end
            end
            for k = 1:numel(this.countriesNames)
                y = this.tProp(:,k);
                b0 = [ones(size(x)) x]\y;
                this.intercept(k) = b0(1);
                this.slope(k) = b0(2);
            end


        end

        function t = getTravellerProportion(this, countryName, year)
            i = find(strcmp(countryName,this.countriesNames));
            if isempty(i)
                t = 1e-6;
            else
                t = max(1e-6, this.intercept(i) + this.slope(i) * year);
            end
        end
    end
end