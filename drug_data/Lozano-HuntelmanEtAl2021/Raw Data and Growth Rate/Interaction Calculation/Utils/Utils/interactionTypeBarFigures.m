function [titleplot,allstats_NN,allstats_EN]=interactionTypeBarFigures(NN,EN,suppCutoff,CI,ndrug,NetOrEmergent)

titleplot=[NetOrEmergent,' ', num2str(ndrug), '-Way Interaction'];

[~,~,~,allstats_NN]=bootstrappercentages(NN,suppCutoff,CI);

%delete 3rd and 4th cols corresponding to buff and supp categorizations
allstats_NN(:,3:4)=[]; % left with 'Synergy','Additive','Antagonism'

barNNdata=cell2mat(allstats_NN(3,:)); %percentages
CI_NN_interactions=cell2mat(allstats_NN(4,:)); %CIs 

NNlowerError=barNNdata-CI_NN_interactions(1,:); %lower error
NNupperError=CI_NN_interactions(2,:)-barNNdata; %upper error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,~,~,allstats_EN]=bootstrappercentages(EN,suppCutoff,CI);

%delete 3rd and 4th cols corresponding to buff and supp categorizations
allstats_EN(:,3:4)=[]; % left with 'Synergy','Additive','Antagonism'

barENdata=cell2mat(allstats_EN(3,:)); %percentages
CI_EN_interactions=cell2mat(allstats_EN(4,:)); %CIs 

ENlowerError=barENdata-CI_EN_interactions(1,:); %lower error
ENupperError=CI_EN_interactions(2,:)-barENdata; %upper error

figure; set(gca,'fontsize',24,'FontName','Arial');
s=size(allstats_EN);
barplotData=[barNNdata;barENdata]';

if isequal(barNNdata,barENdata) %meaning that it is 2-drug data as N2=E2
    NNlowerError=ENlowerError;
    NNupperError=ENupperError;
end

errYbar = zeros(s(2),2,2); 
errYbar(:,:,1)= [NNlowerError; ENlowerError]';
errYbar(:,:,2)=[NNupperError; ENupperError]';

h = barwitherr(errYbar, barplotData); hold on;
set(h(1),'FaceColor','w','EdgeColor','k'); %Net Interaction
set(h(2),'FaceColor','k','EdgeColor','w'); %Emergent Interaction

ymax=ceil(max(barNNdata(:))/0.10)*0.10; ymax=1;
ylim([0,ymax]); set(gca, 'ytick', 0:0.2:ymax); set(gca,'ygrid','on');
xlim([0.5,s(2)+0.5]);
ylabel('Proportion of interactions'); 
title(titleplot,'FontWeight','normal')
set(gca, 'XTick', [1:s(2)], 'XTickLabel', {'Synergy','Additivity','Antagonism'}, 'fontsize',20);
fH = gcf;

end
