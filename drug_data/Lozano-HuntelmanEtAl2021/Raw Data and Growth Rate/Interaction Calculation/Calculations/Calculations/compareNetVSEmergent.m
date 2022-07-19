function [NN_and_ENdata,numberOfCombinationsIncludedInStudy,numberOfCombinationsExcludedFromStudy]=...
                        compareNetVSEmergent(NNdata,ENdata)

%labelsdata={'combinations',[ndrug,'-drug combo'],[ndrug,'-drug combo0'],[NetorEmergent,ndrug],'Interaction'};

basedOnConcentration=1;
ix_values=4; ix_interaction=5;
if basedOnConcentration, ix_Names=2;
else ix_Names=3;
end

NNdataNames=NNdata(:,ix_Names); ENdataNames=ENdata(:,ix_Names);
alldataNames=sort(unique([NNdataNames; ENdataNames]));
ndata=length(alldataNames);

NN_and_ENdata=cell(ndata,3); %combonames, NNdata, NNInteraction, ENdata, ENInteraction
NN_and_ENdata(:,1)=alldataNames; 

% fill the entries with default results: Inconclusive and NaN
NN_and_ENdata(:,3)=repmat({'Inconclusive'},ndata,1); NN_and_ENdata(:,5)=repmat({'Inconclusive'},ndata,1); 
NN_and_ENdata(:,2)=num2cell(nan(ndata,1)); NN_and_ENdata(:,4)=num2cell(nan(ndata,1));


for i=1:ndata
    comboName=alldataNames{i,1};
    ix_NN=find(strcmp(NNdataNames,comboName));
    ix_EN=find(strcmp(ENdataNames,comboName));
    
    if ~isempty(ix_NN) %meaning that it is not inconclusive
        NN_and_ENdata{i,2}=NNdata{ix_NN,ix_values};
        NN_and_ENdata{i,3}=NNdata{ix_NN,ix_interaction};
    end
    
    if ~isempty(ix_EN) %meaning that it is not inconclusive
        NN_and_ENdata{i,4}=ENdata{ix_EN,ix_values};
        NN_and_ENdata{i,5}=ENdata{ix_EN,ix_interaction}; 
    end 
    
end

%find cases with non-inconclusive in both NN and EN

bothAreInconclusive=sum(and(strcmp(NN_and_ENdata(:,3),'Inconclusive'),strcmp(NN_and_ENdata(:,5),'Inconclusive')));
atLeastOneIsNotInconclusive=ndata-bothAreInconclusive;

numberOfCombinationsIncludedInStudy=atLeastOneIsNotInconclusive; %( excluding cases with both Emergent and Net is inconclusive )
numberOfCombinationsExcludedFromStudy=bothAreInconclusive;


end