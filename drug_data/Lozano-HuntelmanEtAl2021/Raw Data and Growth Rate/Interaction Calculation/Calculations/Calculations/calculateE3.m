function [E3,E3_0,E3_1,E3_2,E3_3,E3_0_R,E3_1_R,E3_2_R,E3_3_R,E3Rescaled]=calculateE3(w1,w2,w3,w12,w13,w23,w123,...
                                                                            CheckInconclusive, rescaleMethod)
                                                                        
    if nargin<9, rescaleMethod=2; end                                                                     
                     
    lethalDef=0.047;
    suppressionCutOff=1.3;
    
    lowerOrderParts=-w1*w23-w2*w13-w3*w12+2*w1*w2*w3;
    
    E3=w123+lowerOrderParts;
    
    if E3<=0 %lethal synergy
        wref=0;
    else 
        if rescaleMethod==0, wref=min([w1,w2,w3]);
        elseif rescaleMethod==1, wref=min([w12,w23,w13]);
        elseif rescaleMethod==2, wref=min([w1,w2,w3,w12,w23,w13]);
        else wref=min([w1*w23,w2*w13,w3*w12]); % rescaleMethod=3
        end
    end
    
    E3Rescaled=E3/abs(wref+lowerOrderParts);
    
        [E3_0,E3_1,E3_2,E3_3]=deal(E3/abs(-w1*w23-w2*w13-w3*w12+2*w1*w2*w3));

        if E3>0 %antagonism
            E3_0=E3/abs(min([w1,w2,w3])-w1*w23-w2*w13-w3*w12+2*w1*w2*w3);
            E3_1=E3/abs(min([w12,w23,w13])-w1*w23-w2*w13-w3*w12+2*w1*w2*w3);
            E3_2=E3/abs(min([w1,w2,w3,w12,w23,w13])-w1*w23-w2*w13-w3*w12+2*w1*w2*w3);
            E3_3=E3/abs(min([w1*w23,w2*w13,w3*w12])-w1*w23-w2*w13-w3*w12+2*w1*w2*w3);
        end
        
    if CheckInconclusive && min([w12,w13,w23])<lethalDef && w123<lethalDef
            [E3_0_R,E3_1_R,E3_2_R,E3_3_R]=deal('Inconclusive'); 
            [E3_0,E3_1,E3_2,E3_3]=deal(NaN);
            E3Rescaled=NaN;
    else 
        E3_0_R=InteractionResult(E3_0,suppressionCutOff);
        E3_1_R=InteractionResult(E3_1,suppressionCutOff);
        E3_2_R=InteractionResult(E3_2,suppressionCutOff);
        E3_3_R=InteractionResult(E3_3,suppressionCutOff);
    end       
end


function E3_Result=InteractionResult(E3_scaled,suppressionCutOff)
    E3_Result='Additive';
    if E3_scaled<-0.5, E3_Result='Synergy'; end 
    if E3_scaled>0.5, E3_Result='Antagonistic Buffering'; end 
    if E3_scaled>suppressionCutOff, E3_Result='Antagonistic Suppression'; end 
end