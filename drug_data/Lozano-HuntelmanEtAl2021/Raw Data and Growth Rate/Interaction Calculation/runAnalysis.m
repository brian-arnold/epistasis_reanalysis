%% DESCRIPTION
% This script is intended to reproduce results in Tekin et al. 2018. 
% It calculates interaction metrics for all data contained in Data/input
% and saves the output files in Data/output
%this script was modified to account for corrected data sets
%
% MATLAB version: R2015a
% Operating System: OS X El Capitan Version 10.11.6

% References: 
% Tekin, Elif, Cynthia White, Tina Manzhu Kang, Nina Singh, Mauricio Cruz-Loya, Robert Damoiseaux, Van M. Savage, and Pamela J. Yeh. "Prevalence and patterns of higher-order drug interactions in Escherichia coli." npj Systems Biology and Applications 4, no. 1 (2018): 31.

%% ANALYSIS SETUP: CONFIGS
clear all; clc;
AnalysisSetup;

%% analyze each data file in the dataPath

fprintf('>>>>>>>>>>>>>>>>>>>>>>> ANALYZING ALL FILES in \n')
fprintf([dataPath, '\n'])

[N2MediansAll,N3MediansAll,E3MediansAll, N4MediansAll,E4MediansAll,N5All,E5All]=deal({});
nfiles = length(csvFiles);

    for fileNumber=1:nfiles
        comboName=csvFiles(fileNumber).name;
        fileName= strcat([dataPath,comboName]);
                
        fprintf('******** PROCESSING FILE %s/%s********\n',num2str(fileNumber), num2str(nfiles))
        fprintf([fileName, '\n'])
        
        data = importData(fileName);
        data4and5Drug=data; 
        
        %identify experiments within dataset by empty rows
        skip=identifyExperiments(data);
 
        %create a new data file for storing lower-order information
        [dataN2,dataN3,dataE3,dataN4,dataE4,dataN5,dataE5]=deal({});

        data4and5Drug(1,end+1:end+length(metricLabels))=metricLabels;
        allLabels=[expLabels,metricLabels];
        
    % calculate interaction metrics for each experiment
    % pairwise (only N2), triple(N3,E3), quadruple (N4,E4)),quinary (N5,E5)

        for i=1:size(skip)-1
                %shift one more row due to label row
                data4and5Drug( skip(i)+2: skip(i+1)-1,1)= {data{skip(i)+1,1}};
                combinationName=data{skip(i)+1,1}; combinationName(combinationName == ' ') = [];
                rstr=regexp(combinationName,'([A-Z#+]+)','match'); rconc=regexp(combinationName,'([0-9]+)','match');
                combinationName=strjoin(strcat(rstr,rconc),'+');
                
                combinationdata=cell2mat(data4and5Drug( skip(i)+2: skip(i+1)-1,2:end));
   
                %interaction metric calculations
                [medianfitnesses,metrics,N2Info,N3Info,E3Info,N4Info,E4Info,N5Info,E5Info]=...
                    interactionMetric(combinationdata,combinationName,...
                                                checkInconclusives,CheckSingleCon,calculationMethod,rescaleMethod);
                
                % fill median fitnesses
                data4and5Drug(skip(i)+1,2:2^5)=num2cell(medianfitnesses(2:end));
                data4and5Drug(skip(i)+1,end-length(metricLabels)+1:end)=metrics;
                          
                %lower-orders
                dataN2=[dataN2;N2Info];
                dataN3=[dataN3;N3Info]; dataE3=[dataE3;E3Info];
                dataN4=[dataN4;N4Info]; dataE4=[dataE4;E4Info];
                dataN5=[dataN5;N5Info]; dataE5=[dataE5;E5Info];
        end % end for loop for processing all the experiments for the fileNumberth data file 
        
        % concatanate labels with the median data 
        dataMetrics=[allLabels;data4and5Drug(skip(1:end-1)+1,:)];
        cell2csv(publishData_folder,['/',comboName,'_PublishData','.csv'],  dataMetrics, ',');
        
        basedOnConcentration=1;
        
        %Delete the repeated lower-order infos 
        [~,~,dataN2Medians,~,~,~,~]=lowerOrderClassifications(dataN2,basedOnConcentration,'N');
        [~,~,dataN3Medians,~,~,~,~]=lowerOrderClassifications(dataN3,basedOnConcentration,'N');
        [~,~,dataE3Medians,~,~,~,~]=lowerOrderClassifications(dataE3,basedOnConcentration,'E');
        [~,~,dataN4Medians,~,~,~,~]=lowerOrderClassifications(dataN4,basedOnConcentration,'N');
        [~,~,dataE4Medians,~,~,~,~]=lowerOrderClassifications(dataE4,basedOnConcentration,'E');
        
        N2MediansAll=[N2MediansAll;dataN2Medians];
        N3MediansAll=[N3MediansAll;dataN3Medians];
        E3MediansAll=[E3MediansAll;dataE3Medians];
        N4MediansAll=[N4MediansAll;dataN4Medians];
        E4MediansAll=[E4MediansAll;dataE4Medians];
        N5All=[N5All;dataN5];
        E5All=[E5All;dataE5];
        
        fprintf('******** DONE PROCESSING FILE %s ********\n',num2str(fileNumber))   
    end

% save calculations in results_folder
cell2csv(results_folder,['/','N2 Analysis','.csv'],  N2MediansAll, ',');
cell2csv(results_folder,['/','N3 Analysis','.csv'],  N3MediansAll, ',');
cell2csv(results_folder,['/','E3 Analysis','.csv'],  E3MediansAll, ',');
cell2csv(results_folder,['/','N4 Analysis','.csv'],  N4MediansAll, ',');
cell2csv(results_folder,['/','E4 Analysis','.csv'],  E4MediansAll, ',');

%%  GET MEDIAN of INTERACTION CALCULATION ACROSS EACH EXPERIMENT: THIS IS USED IN PAPER

fprintf('>>>>>>>>>>>>>>>>>>>>>>> GETTING MEDIAN of INTERACTION CALCULATIONS \n')

[~,~,newN2MediansAcrossAll,~,~,~,~]=lowerOrderClassifications(N2MediansAll,basedOnConcentration,'N');
[~,~,newN3MediansAcrossAll,~,~,~,~]=lowerOrderClassifications(N3MediansAll,basedOnConcentration,'N');
[~,~,newE3MediansAcrossAll,~,~,~,~]=lowerOrderClassifications(E3MediansAll,basedOnConcentration,'E');
[~,~,newN4MediansAcrossAll,~,~,~,~]=lowerOrderClassifications(N4MediansAll,basedOnConcentration,'N');
[~,~,newE4MediansAcrossAll,~,~,~,~]=lowerOrderClassifications(E4MediansAll,basedOnConcentration,'E');

s=size(E5All);
for i=1:s(1)
    E5All{i,end}=getInteractionType(E5All{i,end-2},cutoffsupp);
    N5All{i,end}=getInteractionType(N5All{i,end-2},cutoffsupp);
end

% save calculations in results_folder
cell2csv(results_folder,['/','N2 Median Analysis.csv'],  newN2MediansAcrossAll, ',');
cell2csv(results_folder,['/','N3 Median Analysis.csv'],  newN3MediansAcrossAll, ',');
cell2csv(results_folder,['/','N4 Median Analysis.csv'],  newN4MediansAcrossAll, ',');
cell2csv(results_folder,['/','E3 Median Analysis.csv'],  newE3MediansAcrossAll, ',');
cell2csv(results_folder,['/','E4 Median Analysis.csv'],  newE4MediansAcrossAll, ',');
cell2csv(results_folder,['/','E5 Analysis.csv'],  E5All, ',');
cell2csv(results_folder,['/','N5 Analysis.csv'],  N5All, ',');  

%% Create variables for the data analysis summary
close all;

% set indices of names, values and interaction within data files of "newN3MediansAcrossAll" etc
basedOnConcentration=1;
ix_values=4; ix_interaction=6;
if basedOnConcentration, ix_names=2;
else ix_names=3;
end

%distribution bin parameters
binInc=0.1; maxbin=3;

%data collection
ALLn=[2,3,4,5;...
      2,3,4,5];

ALLDATA={newN2MediansAcrossAll,newN3MediansAcrossAll,newN4MediansAcrossAll,N5All;...
         newN2MediansAcrossAll,newE3MediansAcrossAll,newE4MediansAcrossAll,E5All};

ALLTYPES={'Net','Net','Net','Net';...
          'Emergent','Emergent','Emergent','Emergent'};
      
sizeALLDATA=size(ALLDATA);


%% find the size of data
sizeofData=cell(2,5); sizeofData{1,1}='N (number of drugs)'; 
sizeofData{2,1}='number of combinations';
for i=1:length(ALLn)
    N=ALLn(1,i); sizeofData{1,i+1}= N;
    netData=ALLDATA{1,i}; size_netData=size(netData);
    emergentData=ALLDATA{2,i}; size_emergentData=size(emergentData);
    sizeofData{2,i+1}= max([size_netData(1),size_emergentData(1)]);
end
cell2csv(results_folder,['/','number of combinations.csv'],  sizeofData, ',');

%% PLOT breakdown Score vs interaction type proportion

fprintf('>>>>>>>>>>>>>>>>>>>>>>> GETTING BREAKDOWN SCORES \n')

BDS_figures_folder=createFolder(figures_folder,['/Breakdown scores (Fig 3b)']);

for typeIX=1 % 1: Net
    type=ALLTYPES{typeIX,1};
    metric_AcrossAll=ALLDATA(typeIX,:);
    net_AcrossAll=ALLDATA(1,:);
    [combinationNames,Metric,BDS,score_breakdown,perc_breakdown]=breakdownScoresByMedians(metric_AcrossAll,net_AcrossAll,type,BDS_figures_folder);
end

%% PLOT distribution of actual data

fprintf('>>>>>>>>>>>>>>>>>>>>>>> PLOTTING DISTRIBUTIONS \n')

figures_folder_distribution=createFolder(figures_folder,'/Distribution (Fig2)');
results_folder_bootstrap=createFolder(results_folder,'/Bootstrap (Fig 2a)');
                          
for i=1: length(ALLn)
    n=ALLn(1,i); 
    plotNNdata=ALLDATA{1,i}; plotNNdata=plotNNdata(:,ix_values);
    plotENdata=ALLDATA{2,i}; plotENdata=plotENdata(:,ix_values);
    
    %distribution figure
    distributionFigure(n,cell2mat(plotNNdata),cell2mat(plotENdata),binInc,maxbin); %title('Median');
    print( sprintf(['%s/',num2str(n),'drug distribution-median_binsize0p2'], figures_folder_distribution),'-dpng','-r300'); close;
    
    %interaction types bar figure
    [titleplot,allstats_NN,allstats_EN]=interactionTypeBarFigures(plotNNdata,plotENdata,cutoffsupp,CI,n,' ');
    print( sprintf(['%s/',titleplot,''], figures_folder_distribution),'-dpng','-r300'); close;
    
    cell2csv(results_folder_bootstrap,['/','stats N',num2str(n),'.csv'],  allstats_NN, ',');
    cell2csv(results_folder_bootstrap,['/','stats E',num2str(n),'.csv'],  allstats_EN, ',');
end


%% Venn diagrams for comparison of Net and Emergent Interactions

fprintf('>>>>>>>>>>>>>>>>>>>>>>> CREATING VENN DIAGRAMS \n')

results_folder_categorization=createFolder(results_folder,['/Venn (Fig 3a)']);
figures_folder_categorization=createFolder(figures_folder,['/S2 Fig']);

close all;
categorization={'Synergy','Antagonism'};
ndrugLabel={'2-drug','3-drug','4-drug','5-drug'};

for categorization_ix=1:length(categorization)
    interactioncategorization=categorization{categorization_ix};
    interactionPercentages_categorization=cell(sizeALLDATA(1)+1,sizeALLDATA(2)+1);
    interactionPercentages_categorization(1,2:end)=ndrugLabel;
    
    for type=1:sizeALLDATA(1) %row of ALLDATA
        interactionPercentages_categorization{type+1,1}=[ALLTYPES{type,1}, ' ',interactioncategorization];
        interactionPercentages_i=zeros(4,4);
        for ndrug_ix=1:sizeALLDATA(2) % col of ALLDATA
            data=ALLDATA{type,ndrug_ix};

            dataValues=data(:,ix_values);
            [n_syn,n_add,n_buff,n_supp]=InteractionStatistics(cell2mat(dataValues), cutoffsupp);
            interactionPercentages_i(:,ndrug_ix)=[n_syn,n_add,n_buff,n_supp]'./sum([n_syn,n_add,n_buff,n_supp]);
        end
        
        perc_plot=ndrugVsInteractionsFigure(interactionPercentages_i,ALLTYPES{type,1},interactioncategorization); hold on;
        set(gca, 'XTick', [1:4],'XTickLabel',ndrugLabel);
        interactionPercentages_categorization(type+1,2:end)=num2cell(perc_plot);
    end 

    %print legends 
    plot(0.2, 0.55,'ko','MarkerSize',20,'MarkerFaceColor','w');
    text(0.2-0.05, 0.55-0.003, cellstr('N   Net Interaction'), 'FontSize', 15,'FontWeight','normal','Color','k');
    
    plot(0.2, 0.50,'ko','MarkerSize',20,'MarkerFaceColor','k');
    text(0.2-0.04, 0.50-0.003, cellstr('E'), 'FontSize', 15,'FontWeight','normal','Color','w');
    text(0.2-0.04, 0.50-0.003, cellstr('     Emergent Interaction'), 'FontSize', 15,'FontWeight','normal','Color','k');
    
    ylim([0,0.6]); box on; hold off;
    
    print( sprintf(['%s/','Net-Emergent vs ndrug_', interactioncategorization], figures_folder_categorization),'-dpng','-r300'); close;
    
    cell2csv(results_folder_categorization,['/',interactioncategorization,' data.csv'],  interactionPercentages_categorization, ',');
end

%% NET VS EMERGENT VENN DIAGRAM

figures_folder_Venn=createFolder(figures_folder,'/Venn (Fig 3a)');
results_folder_Venn=createFolder(results_folder,'/Venn (Fig 3a)');
        
[NetImpliesEmergentndrug,EmergentImpliesNetndrug]=deal(cell(3,sizeALLDATA(2)+1));
NetImpliesEmergentndrug(:,1)={'NetImpliesEmergent';'Synergy';'Antagonism'};
EmergentImpliesNetndrug(:,1)={'EmergentImpliesNet';'Synergy';'Antagonism'};
NetImpliesEmergentndrug(1,2:end)=ndrugLabel; EmergentImpliesNetndrug(1,2:end)=ndrugLabel; 

for ndrug_ix=1:4 %2-,3-,4-,5-drug
    N=ALLn(1,ndrug_ix); % N drugs in combo
    NNdata=ALLDATA{1,ndrug_ix};
    ENdata=ALLDATA{2,ndrug_ix};
    
    [NN_and_ENdata,numberOfCombinationsIncludedInStudy,numberOfCombinationsExcludedFromStudy]=compareNetVSEmergent(NNdata,ENdata);
    [vennstatsSyn,vennstatsAnt,vennstatsSupp,...
    NetImpliesEmergent_Syn,NetImpliesEmergent_Ant,NetImpliesEmergent_Supp,...
    EmergentImpliesNet_Syn,EmergentImpliesNet_Ant,EmergentImpliesNet_Supp]=...
                                            NetvsEmergent(NN_and_ENdata(:,3), NN_and_ENdata(:,5));
    NetImpliesEmergentndrug{2,ndrug_ix+1}=NetImpliesEmergent_Syn;
    NetImpliesEmergentndrug{3,ndrug_ix+1}=NetImpliesEmergent_Ant;
    
    EmergentImpliesNetndrug{2,ndrug_ix+1}=EmergentImpliesNet_Syn;
    EmergentImpliesNetndrug{3,ndrug_ix+1}=EmergentImpliesNet_Ant;
    
    perc_vennstatsSyn=strcat(cellfun(@(c) num2str(c), num2cell(round(100*(vennstatsSyn/numberOfCombinationsIncludedInStudy))),'UniformOutput', false),'%');
    perc_vennstatsAnt=strcat(cellfun(@(c) num2str(c), num2cell(round(100*(vennstatsAnt/numberOfCombinationsIncludedInStudy))),'UniformOutput', false),'%');
    
    synTitle=[num2str(N),'-Way Synergy']; 
    antTitle=[num2str(N),'-Way Antagonism'];
    
    n_vennstatsSyn=strcat('(n=',cellfun(@(c) num2str(c), num2cell(vennstatsSyn),'UniformOutput', false),')');
    n_vennstatsAnt=strcat('(n=',cellfun(@(c) num2str(c), num2cell(vennstatsAnt),'UniformOutput', false),')');
    
    NetvsEmergentVenn(perc_vennstatsSyn,synTitle,n_vennstatsSyn);
    print( sprintf(['%s/','Net vs Emergent perc', synTitle], figures_folder_Venn),'-dpng','-r300'); close;
    
    NetvsEmergentVenn(perc_vennstatsAnt,antTitle,n_vennstatsAnt);
    print( sprintf(['%s/','Net vs Emergent perc', antTitle], figures_folder_Venn),'-dpng','-r300'); close;

end

cell2csv(results_folder_Venn,'/NetImpliesEmergent data.csv',  NetImpliesEmergentndrug, ',');
cell2csv(results_folder_Venn,'/EmergentImpliesNet data.csv',  EmergentImpliesNetndrug, ',');
