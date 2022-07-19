function skip=identifyExperiments(data)

% this function identifies different experiments from the associated data
% file: data. This is determined based on the rows with blank or NaN
% entries.
% skip: indices of rows that separate different experiments.

    s=size(data); n_rows=s(1);
    
    %find NaN entries in data
    nan_ix=cellfun(@(data) any(isnan(data)),data);
    
    %find empty entries in data
    empty_ix=cellfun('isempty',data);
    
    %find NaN or empty entries in data
    nanORempty_ix=or(nan_ix,empty_ix);
    
    %get rows filled with NaN or empty 
    skip=all(nanORempty_ix,2); skip=find(skip);
    
    %add 0 and n_rows+1 to include start and end the data     
    skip=[0;skip;n_rows+1];
    
    %deal with successive rows of nanORempty_ix in data
    skip(diff(skip)==1)=[];
end