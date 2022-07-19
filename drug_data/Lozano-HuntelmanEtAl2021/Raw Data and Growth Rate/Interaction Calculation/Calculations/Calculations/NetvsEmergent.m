function [vennstatsSyn,vennstatsAnt,vennstatsSupp,...
    NetImpliesEmergent_Syn,NetImpliesEmergent_Ant,NetImpliesEmergent_Supp,...
    EmergentImpliesNet_Syn,EmergentImpliesNet_Ant,EmergentImpliesNet_Supp]=...
                                                                           NetvsEmergent(NNResult, ENResult)

types={'Synergy',{'Antagonistic','Antagonistic Buffering','Antagonistic Suppression'},'Antagonistic Suppression'};
vennstats=cell(1,length(types)); %onlyNN, both, onlyEN
[vennNetImpliesEmergent,vennEmergentImpliesNet]=deal(cell(1,length(types)));


for i=1:length(types)
    type=types{i};
    n_NN=sum(ismember(NNResult,type));
    
    n_EN=sum(ismember(ENResult,type));
    
    n_both=sum(and(ismember(NNResult,type),ismember(ENResult,type)));
    
    n_onlyDA=n_NN-n_both;
    n_onlyEN=n_EN-n_both;
    vennstats{i}=[n_onlyDA,n_both, n_onlyEN];
    
    NetImpliesEmergent=(n_both/n_NN)*100; vennNetImpliesEmergent{i}=NetImpliesEmergent;
    EmergentImpliesNet=(n_both/n_EN)*100; vennEmergentImpliesNet{i}=EmergentImpliesNet;
end

vennstatsSyn=vennstats{1}; vennstatsAnt=vennstats{2}; vennstatsSupp=vennstats{3};
[NetImpliesEmergent_Syn,NetImpliesEmergent_Ant,NetImpliesEmergent_Supp]=...
    deal(vennNetImpliesEmergent{1},vennNetImpliesEmergent{2},vennNetImpliesEmergent{3});

[EmergentImpliesNet_Syn,EmergentImpliesNet_Ant,EmergentImpliesNet_Supp]=...
    deal(vennEmergentImpliesNet{1},vennEmergentImpliesNet{2},vennEmergentImpliesNet{3});

end