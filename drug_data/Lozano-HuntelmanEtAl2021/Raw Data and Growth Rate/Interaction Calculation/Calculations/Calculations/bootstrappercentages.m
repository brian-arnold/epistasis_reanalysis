function [n_supp_data,per_supp_data,CI_supp,allstats]=bootstrappercentages(data,suppCutoff,CI)
% bootstrap the results

%interactionType: Synergy, Additive, Antagonistic Buffering, Antagonistic Suppression,

%this can be implemented by constructing a number of resamples
%with replacement, of the observed dataset (and of equal size to the observed dataset).

data=cell2mat(data); n_tot=sum(~isnan(data));

[n_syn_data,n_add_data,n_buff_data,n_supp_data]=InteractionStatistics(data, suppCutoff);
n_ant_data=n_buff_data+n_supp_data;

%n_supp_data=sum(data>=suppCutoff); updated on April 20,2017
per_supp_data=n_supp_data/n_tot;
M=100;
[n_syn_M,n_add_M,n_buff_M,n_supp_M,n_ant_M]=deal(zeros(1,M));
for i=1:M
    datasample=randsample(data,length(data),1);
    %find interactions
    [n_syn,n_add,n_buff,n_supp]=InteractionStatistics(datasample, suppCutoff);
    %n_supp=sum(datasample>=suppCutoff); 
    n_syn_M(i)=n_syn;
    n_add_M(i)=n_add;
    n_buff_M(i)=n_buff;
    n_supp_M(i)=n_supp;
    n_ant_M(i)=n_buff+n_supp;
end

n_syn_M=sort(n_syn_M); n_add_M=sort(n_add_M); 
n_buff_M=sort(n_buff_M); n_supp_M=sort(n_supp_M); n_ant_M=sort(n_ant_M);

alpha=(100-CI)/2;
Lower_perc_ix=floor(alpha*M/100);
Upper_per_ix=ceil((100-alpha)*M/100);

%n_syn_M

CI_syn=[n_syn_M(Lower_perc_ix),n_syn_M(Upper_per_ix)]; CI_syn=CI_syn/n_tot;
CI_add=[n_add_M(Lower_perc_ix),n_add_M(Upper_per_ix)]; CI_add=CI_add/n_tot;
CI_buff=[n_buff_M(Lower_perc_ix),n_buff_M(Upper_per_ix)]; CI_buff=CI_buff/n_tot;
CI_supp=[n_supp_M(Lower_perc_ix),n_supp_M(Upper_per_ix)]; CI_supp=CI_supp/n_tot;
CI_ant=[n_ant_M(Lower_perc_ix),n_ant_M(Upper_per_ix)]; CI_ant=CI_ant/n_tot;


allstats={'Synergy','Additive','Antagonistic Buffering','Antagonistic Suppression','Antagonism';...
    n_syn_data,n_add_data,n_buff_data,n_supp_data,n_ant_data;...
    n_syn_data/n_tot,n_add_data/n_tot,n_buff_data/n_tot,n_supp_data/n_tot,n_ant_data/n_tot;...
    CI_syn',CI_add',CI_buff',CI_supp',CI_ant'};

end