function [n_syn,n_add,n_buff,n_supp]=InteractionStatistics(metric_results, suppCutoff)

%synergy between -1,-0.5
syn_metrics=metric_results(metric_results<-0.5);
n_syn=length(syn_metrics);

%additive between -0.5,0.5
add_metrics=metric_results(and(metric_results<0.5,metric_results>-0.5));
n_add=length(add_metrics);

%additive between 0.5,suppCutoff
buff_metrics=metric_results(and(metric_results<suppCutoff,metric_results>0.5));
n_buff=length(buff_metrics);

%suppression between suppCutoff,inf
supp_metrics=metric_results(metric_results>=suppCutoff);
n_supp=length(supp_metrics);

end