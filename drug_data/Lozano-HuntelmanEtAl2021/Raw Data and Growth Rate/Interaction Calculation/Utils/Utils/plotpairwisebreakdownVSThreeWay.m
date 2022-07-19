function plotpairwisebreakdownVSThreeWay(score_breakdown,perc_breakdown)
    types={'Synergy','Additive','Antagonistic'};
    colors={'r','k',[0 153 51]./255};
    
    for i=1:length(types)
        if i~=2 %Skip additive curve
            plot(score_breakdown,perc_breakdown(:,i),'.-','Color',colors{i},'MarkerSize',20); hold on;
            ylabel('Proportion of interaction type','FontWeight','normal'); xlabel('breakdown score');
            
            xmax=score_breakdown(end); 
            xlim([-xmax,xmax]);
            set(gca,'XTick',[-xmax,0,xmax])
            axis square;
        end
    end 
    
    set(gca,'fontsize',20, 'FontName','Arial');
end