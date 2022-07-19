function [per_supp_MedianEpsilons,per_supp_AllEpsilons,...
    newdataMedians,newdata,sortednewdata,uniquenewdata,values]=lowerOrderClassifications(dataLower,basedOnConcentration,NetorEmergent)

    uniquenewdata=[];

    ix_values=4; ix_interaction=6; ix_fitness=5;
    if basedOnConcentration, ix_Names=2;
    else ix_Names=3;
    end
    
    %sort lower-order data: same combos are ordered successively.
    % NOTE: combination names in dataLower 
    newdata=sortPairwiseData(dataLower,basedOnConcentration);
    sortednewdata=newdata;
    
    %delete duplicate experiments
    values=cell2mat(newdata(:,ix_values));
    
    newdata=newdata(~isnan(values),:);
    
    values=cell2mat(newdata(:,ix_values));
    doubleNames=newdata(:,ix_Names);
    num_dig=5; values = round(values.*(10^num_dig))./(10^num_dig); %done to get over with the roundness issues!
    newdata(:,ix_values)=num2cell(values);
    
    ndrug=num2str(length(strsplit(doubleNames{1},'+')));
    
    
    cutoffsupp=1.3;
    if ndrug==2, cutoffsupp=1.15; end 
    
    %find indices of unique pairwise combinations 
    [name_unique,ix_unique_name]=unique(doubleNames);
    
    %pairwise comb names
    newdataMedians=newdata(ix_unique_name,:);
    s=size(newdata); n_rows=s(1);
    
    ix_unique_name=[ix_unique_name;n_rows+1];
    for i=1:length(ix_unique_name)-1
        newdata{ix_unique_name(i),4};
        
        DAs=unique(cell2mat(newdata(ix_unique_name(i):ix_unique_name(i+1)-1,ix_values)));
        medianDA=nanmedian(DAs);
        newdataMedians{i,ix_values}=medianDA;
        
        Ws=unique(cell2mat(newdata(ix_unique_name(i):ix_unique_name(i+1)-1,ix_fitness)));
        medianW=nanmedian(Ws);
        newdataMedians{i,ix_fitness}=medianW;
        
        newdataMedians{i,ix_interaction}=getInteractionType(medianDA,cutoffsupp);
    end
    
    [per_supp_MedianEpsilons,~]=getInteractionClassificationNumbers(newdataMedians(:,ix_interaction));
    [per_supp_AllEpsilons,~]=getInteractionClassificationNumbers(newdata(:,ix_interaction));
    
    %add labels on top
    labelsPairwise={'combinations',[ndrug,'-drug combo'],[ndrug,'-drug combo0'],[NetorEmergent,ndrug], 'Fitness','Interaction'};
    newdata=[labelsPairwise;newdata]; 
    
end