function perc_plot=ndrugVsInteractionsFigure(interactionPercentages,NetOrEmergent, interactionType)

    perc_syn=interactionPercentages(1,:); 
    perc_add=interactionPercentages(2,:); 
    perc_buff=interactionPercentages(3,:);
    perc_supp=interactionPercentages(4,:);
    perc_ant=perc_buff+perc_supp;
    
    YLabel=interactionType;
    
    if strcmp(interactionType,'Synergy')
        perc_plot=perc_syn;
    elseif strcmp(interactionType,'Additivity')
        perc_plot=perc_add;
    elseif strcmp(interactionType,'Antagonism')
        perc_plot=perc_ant;
    else 
        perc_plot=perc_supp;
    end
    
    plotColor='w'; textColor='k';
    if strcmp(NetOrEmergent,'Emergent'), plotColor='k'; textColor='w'; end
    
    s=size(interactionPercentages);
    xvalues=1:s(2); 
    
    %figure; 
    set(gca,'fontsize',20,'FontName','Arial'); hold on;
    
    plot(xvalues, perc_plot,'ko','MarkerSize',20,'MarkerFaceColor',plotColor,'LineStyle','-');
    
    for i=1:s(2)
        text(xvalues(i)-0.04, perc_plot(i)-0.001, cellstr(NetOrEmergent(1)), 'FontSize', 15,'FontWeight','normal','Color',textColor);
    end
   
    ymax=1; ylim([0,ymax]);
    set(gca, 'ytick', 0:0.1:ymax); 
    set(gca,'ygrid','on');  ylabel(['Proportion of ', YLabel]); 
    xlim([0,s(2)+0.2]); set(gca, 'XTick', [1:s(2)])
   
    b = bar(0,'EdgeColor','none');
end