function [percentage_supp,n_interactions]=getInteractionClassificationNumbers(metric_results)

    n_supp=sum(strcmp(metric_results,'Antagonistic Suppression'));
    n_buff=sum(strcmp(metric_results,'Antagonistic'))+sum(strcmp(metric_results,'Antagonistic Buffering'));
    n_syn=sum(strcmp(metric_results,'Synergy'));
    n_add=sum(strcmp(metric_results,'Additive'));
    n_inc=sum(strcmp(metric_results,'Inconclusive'));
    n_tot=length(metric_results);
    n_interactions={'n_tot', 'n_inc','n_syn','n_add','n_buff','n_supp';...
    n_tot,n_inc,n_syn,n_add,n_buff,n_supp};   
    percentage_supp=n_supp/(n_tot-n_inc);

end