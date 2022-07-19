function interaction=getInteractionType(metric_rescaled,cutoffsupp)
    interaction='Additive';
    if isnan(metric_rescaled), interaction='Inconclusive'; end
    if metric_rescaled<-0.5, interaction='Synergy';  end
    if metric_rescaled>0.5 && metric_rescaled<cutoffsupp, interaction='Antagonistic'; end
    if metric_rescaled>cutoffsupp, interaction='Antagonistic Suppression'; end
end