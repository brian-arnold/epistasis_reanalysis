function [N3,N3_Scaled,N3_Result,N3_fit]=calculateN3(wX,wY,wZ,wXY,wXZ,wYZ,wXYZ,CheckInconclusive)
    N3=wXYZ-wX*wY*wZ; 
    N3_Scaled=N3/(wX*wY*wZ); if N3>0, N3_Scaled=N3/abs(min([wX,wY,wZ])-wX*wY*wZ); end
    N3_fit = wXYZ;
    suppressionCutOff=1.3;
    N3_Result='Additive';
    if N3_Scaled<-0.5, N3_Result='Synergy'; end 
    if N3_Scaled>0.5, N3_Result='Antagonistic Buffering'; end
    if N3_Scaled>suppressionCutOff, N3_Result='Antagonistic Suppression'; end 
    
    
    lethalDef=0.04; lethalDef=0.047; WTDef=0.9;
    
    if nargin == 7   % if the number of inputs equals 2
         CheckInconclusive = 0; % then make the third value, z, equal to my default value, 5.
    end
    
    if CheckInconclusive && min([wX,wY,wZ])<lethalDef && wXYZ<lethalDef
           N3_Scaled=NaN; N3_Result='Inconclusive'; 
            N3_fit = NaN;
    end    
    
    
end