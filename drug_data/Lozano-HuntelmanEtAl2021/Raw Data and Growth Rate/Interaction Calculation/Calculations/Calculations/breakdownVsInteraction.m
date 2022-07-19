function [combinationNames,Metric,BDS,score_breakdown,perc_breakdown]=breakdownVsInteraction(MetricvsBreakdownScore,N)
    
    combinationNames=MetricvsBreakdownScore(:,2);
    MetricvsBreakdownScore=cell2mat(MetricvsBreakdownScore(:,end-1:end));
    Metric=MetricvsBreakdownScore(:,1);
    BDS=MetricvsBreakdownScore(:,2);
    
    [~,uniqueIndicesMetric]=unique(Metric); [~,uniqueIndicesBDS]=unique(BDS);
    [~,uniqueIndicesNames]=unique(combinationNames);
    uniqueIndices=intersect(uniqueIndicesMetric,uniqueIndicesNames);
    
    %NOTE: because MetricvsBreakdownScore is medians across all, these are
    %already unique values
    if length(uniqueIndicesNames)~=length(combinationNames)
        length(uniqueIndicesNames)
        length(combinationNames)
        fprintf('something is wrong, check breakdownVsInteraction.m file');
    end
    uniqueIndices=uniqueIndicesNames;
    combinationNames=combinationNames(uniqueIndices);
    Metric=Metric(uniqueIndices);
    BDS=BDS(uniqueIndices);


    maxScore=0;
    for i=N-1:-1:2
        maxScore=maxScore+nchoosek(N,i);
    end 

    if max(BDS)> maxScore || min(BDS)< -maxScore
        error('max break down is not correct, check breakdownVsInteraction.m')
    end
    
    score_breakdown=[-maxScore:1:maxScore]';
    [AllVector,SynergyVector, AntagonismVector,AdditiveVector]=deal(zeros(size(score_breakdown)));

    ix=1;
    for score=-maxScore:maxScore
        AllCount=sum(BDS==score); AllVector(ix)=AllCount;
        SynergyCount=sum(Metric(BDS==score)<-0.5); SynergyVector(ix)=SynergyCount;
        AntagonismCount=sum(Metric(BDS==score)>0.5); AntagonismVector(ix)=AntagonismCount;
        AdditiveCount=AllCount-SynergyCount-AntagonismCount; AdditiveVector(ix)=AdditiveCount;
        ix=ix+1;
    end
    SynergyProportion=SynergyVector./AllVector;
    AntagonismProportion=AntagonismVector./AllVector;
    AdditiveProportion=AdditiveVector./AllVector;
    
    perc_breakdown=[SynergyProportion,AdditiveProportion,AntagonismProportion];

end