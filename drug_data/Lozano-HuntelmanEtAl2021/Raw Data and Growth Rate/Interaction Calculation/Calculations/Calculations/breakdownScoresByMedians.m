function [combinationNames,Metric,BDS,score_breakdown,perc_breakdown]=breakdownScoresByMedians(metric_AcrossAll,net_AcrossAll,type, folder)
%for any N drug combo, find all the lower-order combinations, and create
%breakdown scores
clc;
    for N=3:5
        %N: number of drugs > 3
        metric_NDrugInfo=metric_AcrossAll{N-1};
        l=length(metric_NDrugInfo);
        metric_NDrug_BDS=zeros(l,1);

        for i=1:l %loop over all Ndrug combos
            NdrugCombo=metric_NDrugInfo{i,2};
            drugs=strsplit(NdrugCombo,'+');
            BDS=0;
            for j=N-1:-1:2 %loop over all lower-order combos of NdrugCombo
                metric_jDrugInfo=net_AcrossAll{j-1};
                metric_jDrugComboNames=metric_jDrugInfo(:,2);
                metric_jDrugComboMetrics=metric_jDrugInfo(:,end-1);

                jdrugcombinations=nchoosek(drugs,j);
                for ix=1:nchoosek(N,j)
                    jdrugComboName=strjoin(sort(jdrugcombinations(ix,:)),'+');
                    indexjDdrugCombo=find(strcmp(metric_jDrugComboNames,jdrugComboName));

                    if  ~isempty(indexjDdrugCombo) %~isempty(jDrugComboMetric) &&
                        jDrugComboMetric=metric_jDrugComboMetrics{indexjDdrugCombo};

                        if jDrugComboMetric<-0.5, BDS=BDS-1;
                        elseif jDrugComboMetric>0.5, BDS=BDS+1;
                        end 

                    else 
                        BDS=NaN;
                    end
                end
            end
            
            metric_NDrug_BDS(i)=BDS;

        end

        metricvsBDS_N=[metric_NDrugInfo(:,1:end-1),num2cell(metric_NDrug_BDS)];

        figure;
        [combinationNames,Metric,BDS,score_breakdown,perc_breakdown]=breakdownVsInteraction(metricvsBDS_N,N);
        plotpairwisebreakdownVSThreeWay(score_breakdown,perc_breakdown); 
        title([type,' ', num2str(N), '-drug']);
        print( sprintf(['%s/',type,' ', num2str(N), '-drug',' vs BDS'], folder),'-dpng','-r300'); close;
        cell2csv(folder,['/', type, num2str(N),' vs BDS_','.csv'],[combinationNames,num2cell(Metric),num2cell(BDS)],',');
    
    end 

end

