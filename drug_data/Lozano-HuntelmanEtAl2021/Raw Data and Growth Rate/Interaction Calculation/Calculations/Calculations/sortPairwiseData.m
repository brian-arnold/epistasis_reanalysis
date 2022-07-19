function newdataLower=sortPairwiseData(dataLower,basedOnConcentration)
    %sort based on pairwise: same pairs are ordered successively!
    newdataLower={};
    if basedOnConcentration, lowerOrderComboNames=dataLower(:,2);
    else lowerOrderComboNames=dataLower(:,3);
    end

    allcombos=unique(lowerOrderComboNames); %get unique combo names
    for i=1: length(allcombos)
        %identify experiments of the combo: allexps{i}
        exp_i=strcmp(lowerOrderComboNames,allcombos{i});
        nextLowerOrderCombo=dataLower(exp_i,:);
        
        %as well as based on the values???
        [~, idx] = sort([nextLowerOrderCombo{:,end-1}], 'ascend'); 
        nextLowerOrderCombo=nextLowerOrderCombo(idx,:);
       
        % add a new line
        newdataLower=[newdataLower;nextLowerOrderCombo];
    end
end
