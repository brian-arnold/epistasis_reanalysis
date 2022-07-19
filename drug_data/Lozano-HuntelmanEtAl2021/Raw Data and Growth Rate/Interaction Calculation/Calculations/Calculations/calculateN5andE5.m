function [N5rescaled,E5Rescaled,N5_fit]=...
                                calculateN5andE5(w1,w2,w3,w4,w5,...
                                w12,w13,w14,w15,w23,w24,w25,w34,w35,w45,...
                                w123,w124,w125,w134,w135,w145,w234,w235,w245,w345,...
                                w1234,w1235,w1245,w1345,w2345,...
                                w12345,...
                                CheckInconclusive, method,rescaleMethod)
%% Constant parameters

lethalDef=0.047; WTDef=0.9;

%% NET INTERACTION
if CheckInconclusive && min([w1,w2,w3,w4,w5])<lethalDef && w12345<lethalDef   
    N5rescaled=NaN;
     N5_fit = NaN;  
else 
    N5=w12345-w1*w2*w3*w4*w5;
    N5rescaled=N5/(w1*w2*w3*w4*w5);

    if N5>0, N5rescaled=N5/(abs(min([w1,w2,w3,w4,w5])-w1*w2*w3*w4*w5)); end 
    N5_fit = w12345;
end

%% EMERGENT INTERACTION

if CheckInconclusive && min([w1234,w1235,w1245,w1345,w2345])<lethalDef && w12345<lethalDef
   E5Rescaled=NaN;
   
else
    
    if strcmp(method,'Ursell')
        lowerOrderParts=-(w1*w2345+w2*w1345+w3*w1245+w4*w1235+w5*w1234)...
                 -(w12*w345+w13*w245+w14*w235+w15*w234+w23*w145+w24*w135+w25*w134+w34*w125+w35*w124+w45*w123)...
                 +2*(w4*w5*w123+w3*w5*w124+w3*w4*w125+w2*w5*w134+w2*w4*w135+w2*w3*w145+w1*w5*w234+w1*w4*w235+w1*w3*w245+w1*w2*w345)...
                 +2*(w12*w34*w5+w12*w35*w4+w12*w45*w3+w13*w24*w5+w13*w25*w4+w13*w45*w2+w14*w23*w5+w14*w25*w3+w14*w35*w2+w15*w23*w4+w15*w24*w3+w15*w34*w2+w23*w45*w1+w24*w35*w1+w25*w34*w1)...
                 -6*(w1*w2*w3*w45+w1*w2*w4*w35+w1*w2*w5*w34+w1*w3*w4*w25+w1*w3*w5*w24+w1*w4*w5*w23+w2*w3*w4*w15+w2*w3*w5*w14+w2*w4*w5*w13+w3*w4*w5*w12)...
                 +24*w1*w2*w3*w4*w5;
    elseif strcmp(method,'Covariance')   
        lowerOrderParts=-w1*w2345-w2*w1345-w3*w1245-w4*w1235-w5*w1234...
             +w4*w5*w123+w3*w5*w124+w3*w4*w125+w2*w5*w134+w2*w4*w135+w2*w3*w145+w1*w5*w234+w1*w4*w235+w1*w3*w245+w1*w2*w345...
             -w1*w2*w3*w45-w1*w2*w4*w35-w1*w2*w5*w34-w1*w3*w4*w25-w1*w3*w5*w24-w1*w4*w5*w23-w2*w3*w4*w15-w2*w3*w5*w14-w2*w4*w5*w13-w3*w4*w5*w12...
             +4*w1*w2*w3*w4*w5;   

    end 
     
    E5=w12345+lowerOrderParts;
    
    if E5<=0 % rescale relative to lethal synergy
        wref=0;
    else % rescale relative to antagonistic buffering
        
        if rescaleMethod==2 % min of lower-order measurements
            wref=min([w1,w2,w3,w4,w5,...
                w12,w13,w14,w15,w23,w24,w25,w34,w35,w45,...
                w123,w124,w125,w134,w135,w145,w234,w235,w245,w345,...
                w1234,w1235,w1245,w1345,w2345]);
        
        elseif rescaleMethod==3  % min of lower-order effects 
            
            if strcmp(method,'Ursell') 
                wref=min([w1*w2345,w2*w1345,w3*w1245,w4*w1235,w5*w1234,...
                     w12*w345,w13*w245,w14*w235,w15*w234,w23*w145,w24*w135,w25*w134,w34*w125,w35*w124,w45*w123,...
                     w4*w5*w123,w3*w5*w124,w3*w4*w125,w2*w5*w134,w2*w4*w135,w2*w3*w145,w1*w5*w234,w1*w4*w235,w1*w3*w245,w1*w2*w345,...
                     w12*w34*w5,w12*w35*w4,w12*w45*w3,w13*w24*w5,w13*w25*w4,w13*w45*w2,w14*w23*w5,w14*w25*w3,w14*w35*w2,w15*w23*w4,w15*w24*w3,w15*w34*w2,w23*w45*w1,w24*w35*w1,w25*w34*w1,...
                     w1*w2*w3*w45,w1*w2*w4*w35,w1*w2*w5*w34,w1*w3*w4*w25,w1*w3*w5*w24,w1*w4*w5*w23,w2*w3*w4*w15,w2*w3*w5*w14,w2*w4*w5*w13,w3*w4*w5*w12]);
            elseif strcmp(method,'Covariance')
                wref=min([w1*w2345,w2*w1345,w3*w1245,w4*w1235,w5*w1234,...
                     w4*w5*w123,w3*w5*w124,w3*w4*w125,w2*w5*w134,w2*w4*w135,w2*w3*w145,w1*w5*w234,w1*w4*w235,w1*w3*w245,w1*w2*w345,...
                     w1*w2*w3*w45,w1*w2*w4*w35,w1*w2*w5*w34,w1*w3*w4*w25,w1*w3*w5*w24,w1*w4*w5*w23,w2*w3*w4*w15,w2*w3*w5*w14,w2*w4*w5*w13,w3*w4*w5*w12]);
            end
            
        end
    end 
    
    E5Rescaled=E5/(abs(wref+lowerOrderParts));
end  


end