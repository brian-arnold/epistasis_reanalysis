function [N2_Scaled,N2_Result,N2_fit]=calculateN2(w1,w2,w12)

    lowerOrderParts=-w1*w2;
    N2=w12+lowerOrderParts;
        
    if N2<=0
        wref=0;
    else %N2>0
        wref=min([w1,w2]);
    end
    
    N2_Scaled=N2/abs(wref+lowerOrderParts); 
    N2_fit = w12;
    
    cutoffsupp=1.15;
    N2_Result=getInteractionType(N2_Scaled,cutoffsupp);
    
    WTDef=0.9;
    if max(w1,w2)>WTDef && w12>WTDef, N2_Scaled=NaN; N2_Result='Inconclusive'; N2_fit = NaN; end 
   
end