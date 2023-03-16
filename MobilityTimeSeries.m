classdef MobilityTimeSeries
    properties
        mobilityIndicators
        years
    end
    methods
        function this = MobilityTimeSeries(mobilityIndicators)
            this.mobilityIndicators = mobilityIndicators;
            clear mobilityIndicators
            for i = 1:numel(this.mobilityIndicators)
                disp(['Storing sMatrix ' num2str(i) ' out of ' num2str(numel(this.mobilityIndicators))])
                this.mobilityIndicators(i) = this.mobilityIndicators(i).store_sMatrix;
                this.years(i) = this.mobilityIndicators(i).year;
            end
        end
    end
    
    methods 
        function perLocationAverageImportsTimeSeries(this, fnum)
            figure(fnum)
            clf(fnum)
            for i = 1:numel(this.mobilityIndicators)
                [~,imports] = perLocationFlowRate(this.mobilityIndicators(i));
%                 imp(:,i) = imports;
                m(i) = mean(imports);
                mn(i) = min(imports);
                mx(i) = max(imports);
                cv(i) = m(i) / std(imports);
                
            end
            x = this.years;
            
            
            y1 = mn;%m - cv/2;
            y2 = mx;%m + cv/2;
%             plot(x', [y1' y2'] ,'-')
            
            p =patch([x fliplr(x)], [y1 fliplr(y2)], 'c');
            set(p,'EdgeColor','none')
            hold on
            plot(x, m,'-k')
            hold off
        end
        function totalImportsTimeSeries(this, fnum)
            figure(fnum)
            clf(fnum)
            for i = 1:numel(this.mobilityIndicators)
                mi = this.mobilityIndicators(i);
                
                [~,imports_CasesWeight] = mi.perLocationFlowRate(mi.travellerProportion);
                infP = mi.dataTable.infectedProportion';
                infP(mi.travellerProportion == 0) = 1;
                [~,imports_NoCasesWeight] = mi.perLocationFlowRate(mi.travellerProportion ./ infP);

%                 imp(:,i) = imports;
                t(i) = sum(imports_CasesWeight / 12);
                tLTS(i) = sum((imports_CasesWeight / 12) .* mi.dataTable(mi.dataTable.logicalImportLocations,:).LTS);
                t2(i) = sum(imports_CasesWeight);
                mLTS(i) = mean( mi.dataTable(mi.dataTable.logicalImportLocations,:).LTS );
                a(i) = sum(mi.dataTable.logicalImportLocations);

                t_(i) = sum(imports_NoCasesWeight / 12);
                tLTS_(i) = sum((imports_NoCasesWeight / 12) .* mi.dataTable(mi.dataTable.logicalImportLocations,:).LTS);
                t2_(i) = sum(imports_NoCasesWeight);
                mLTS_(i) = mean( mi.dataTable(mi.dataTable.logicalImportLocations,:).LTS );
                a_(i) = sum(mi.dataTable.logicalImportLocations);

            end

            x = this.years;
            ax(1) = subplot(3,2,1);
            plot(x, tLTS,'-m')
            hold on
            plot(x, this.lin(x,tLTS),'-k')
            ylabel('Yearly imports (Monthly imports * LTS)')

            ax(2) = subplot(3,2,3);
            plot(x, mLTS,'-m')
            hold on
            plot(x, this.lin(x,mLTS),'-k')
            ylabel('mean LTS')
            
            ax(3) = subplot(3,2,5);
            plot(x, a,'-m')
            hold on
            plot(x, this.lin(x,a),'-k')
            ylabel('Number of NUTS-regions with LTS > 1 month')

            ax(4) = subplot(3,2,2);
            plot(x, tLTS_,'-m')
            hold on
            plot(x, this.lin(x,tLTS_),'-k')
            ylabel('Yearly number of visitors (Monthly number of visitors * LTS')

            ax(5) = subplot(3,2,4);
            plot(x, mLTS_,'-m')
            hold on
            plot(x, this.lin(x,mLTS_),'-k')
            ylabel('mean LTS')
            
            ax(6) = subplot(3,2,6);
            plot(x, a_,'-m')
            hold on
            plot(x, this.lin(x,a_),'-k')
            ylabel('Number of NUTS-regions with LTS > 1 month')
            
            set(gca,'YScale','linear')
            xlim(ax,[min(this.years) max(this.years)])
            hold off
        end

        function [X,Y] = mobilityTrendFigureXY(this, aggregationOperator)
            arguments
                this MobilityTimeSeries
                aggregationOperator string
            end
            X = this.years;
            for i = 1:numel(this.mobilityIndicators)
                mi = this.mobilityIndicators(i);
                
                [~,imports_CasesWeight] = mi.perLocationFlowRate(mi.travellerProportion);
%                 infP = mi.dataTable.infectedProportion';
%                 infP(mi.travellerProportion == 0) = 1;
%                 [~,imports_NoCasesWeight] = mi.perLocationFlowRate(mi.travellerProportion ./ infP);
                if strcmpi(aggregationOperator,'sum')
                    Y(i) = sum((imports_CasesWeight / 12) .* mi.dataTable(mi.dataTable.logicalImportLocations,:).LTS);
                elseif strcmpi(aggregationOperator,'mean')
                    Y(i) = mean((imports_CasesWeight / 12) .* mi.dataTable(mi.dataTable.logicalImportLocations,:).LTS);
                else
                    error('Use either "sum" or "mean" as aggregationOperator.')
                end
            end
        end
    end

    methods
        function Y = lin(this,x,y)
            arguments
                this MobilityTimeSeries
                x (:,1) double
                y (:,1) double
            end
            X = [ones(size(x)) x];
            b = X\y;
            Y = X*b;
        end

    end
end