function strNew=deleteConc(str)
    r=regexp(str,'([A-Z#+]+)','match');
    strNew=strjoin(r,'');
end