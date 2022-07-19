function [N4rescaled,E4Rescaled, N4_fit]=calculateN4andE4(w1,w2,w3,w4,...
                                    w12,w13,w14,w23,w24,w34,...
                                    w123,w124,w134,w234,...
                                    w1234,...
                                    CheckInconclusive,method,rescaleMethod)
%% Constant parameters

lethalDef=0.047;
WTDef=0.9;

%% NET INTERACTION
    
if CheckInconclusive && min([w1,w2,w3,w4])<lethalDef && w1234<lethalDef
    N4rescaled=NaN;
    N4_fit = NaN;
else
    N4=w1234-w1*w2*w3*w4;
    N4rescaled=N4/(w1*w2*w3*w4);
    if N4>0, N4rescaled=N4/(abs(min([w1,w2,w3,w4])-w1*w2*w3*w4)); end 
    N4_fit = w1234;
end    


%% EMERGENT INTERACTION
if CheckInconclusive && min([w123,w124,w134,w234])<lethalDef && w1234<lethalDef
  
    E4Rescaled=NaN;
       
else 
    if strcmp(method,'Ursell')
        lowerOrderParts=-(w1*w234+w2*w134+w3*w124+w4*w123)...
                    -(w12*w34+w13*w24+w14*w23)...
                    +2*(w1*w2*w34+w1*w3*w24+w1*w4*w23+w2*w3*w14+w2*w4*w13+w3*w4*w12)...
                    -6*w1*w2*w3*w4;
    elseif strcmp(method,'Covariance')
        lowerOrderParts=-w1*w234-w2*w134-w3*w124-w4*w123...
                    +w1*w2*w34+w1*w3*w24+w1*w4*w23+w2*w3*w14+w2*w4*w13+w3*w4*w12...
                    -3*w1*w2*w3*w4;         
    elseif strcmp(method,'Isserlis')
        lowerOrderParts=-(w12*w34+w13*w24+w14*w23)+2*w1*w2*w3*w4;
    end

    E4=w1234+lowerOrderParts;

    if E4<=0
        wref=0;
    else
        if rescaleMethod==2 
            
            if strcmp(method,'Ursell') || strcmp(method,'Covariance') 
                wref=min([w1,w2,w3,w4,w12,w13,w14,w23,w24,w34,w123,w124,w134,w234]);
            elseif strcmp(method,'Isserlis')
                wref=min([w12,w13,w14,w23,w24,w34]);
            end
            
        elseif rescaleMethod==3
            
            if strcmp(method,'Ursell')
                wref=min([w1*w234,w2*w134,w3*w124,w4*w123,...
                        w12*w34,w13*w24,w14*w23,...
                        w1*w2*w34,w1*w3*w24,w1*w4*w23,w2*w3*w14,w2*w4*w13,w3*w4*w12]);
            elseif strcmp(method,'Covariance')
                wref=min([w1*w234,w2*w134,w3*w124,w4*w123,...
                        w1*w2*w34,w1*w3*w24,w1*w4*w23,w2*w3*w14,w2*w4*w13,w3*w4*w12]);
            elseif strcmp(method,'Isserlis')
                wref=min([w12*w34,w13*w24,w14*w23]);
            end 
            
        end
    end 
    
    E4Rescaled=E4/(abs(wref+lowerOrderParts));
end



end