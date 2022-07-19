function [med_w,metrics,N2Info,N3Info,E3Info,N4Info,E4Info,N5Info,E5Info]=...
                            interactionMetric(combinationdata,combinationName,...
                                                        CheckInconclusive,CheckSingleCon,...
                                                        method,rescaleMethod)
    alpha=0.05; 
    singleNames=strsplit(combinationName,'+');
    nsingles=length(singleNames);
    s=size(combinationdata); combinationdata=[100*ones(s(1),1),combinationdata];
    
    combinationdata (combinationdata==-111)=NaN; % this for cleaning extra -111s in the data file that were place holders
    med_w=nanmedian(combinationdata);
    
    med_w(med_w>100)=100; med_w(med_w<0)=0; 
    
    med_w=med_w/100; %convert to percentage values, hence between 0 and 1.
   
    [w0,w1,w2,w3,w4,w5,...
     w12,w13,w14,w15,w23,w24,w25,w34,w35,w45,...
     w123,w124,w125,w134,w135,w145,w234,w235,w245,w345,...
     w1234,w1235,w1245,w1345,w2345,...
     w12345]=deal(med_w(1),med_w(2),med_w(3), med_w(4),med_w(5),...
                  med_w(6),med_w(7),med_w(8),med_w(9),med_w(10),med_w(11),med_w(12),med_w(13),med_w(14),med_w(15),...
                  med_w(16),med_w(17),med_w(18),med_w(19),med_w(20),med_w(21),med_w(22),med_w(23),med_w(24),med_w(25),...
                  med_w(26),med_w(27),med_w(28),med_w(29),med_w(30),...
                  med_w(31),med_w(32));
 
%% 2-way interactions             
    [DA_12,~,w_12]=calculateN2(w1,w2,w12);
    [DA_13,~,w_13]=calculateN2(w1,w3,w13);
    [DA_14,~,w_14]=calculateN2(w1,w4,w14);
    [DA_15,~,w_15]=calculateN2(w1,w5,w15);
    [DA_23,~,w_23]=calculateN2(w2,w3,w23);
    [DA_24,~,w_24]=calculateN2(w2,w4,w24);
    [DA_25,~,w_25]=calculateN2(w2,w5,w25);
    [DA_34,~,w_34]=calculateN2(w3,w4,w34);
    [DA_35,~,w_35]=calculateN2(w3,w5,w35);
    [DA_45,~,w_45]=calculateN2(w4,w5,w45);
    
    DA2=[DA_12,DA_13,DA_14,DA_15,DA_23,DA_24,DA_25,DA_34,DA_35,DA_45];
    W2=[w_12,w_13,w_14, w_15,w_23,w_24,w_25,w_34,w_35,w_45];
    
    if CheckSingleCon && max([w1,w2,w3,w4,w5])>0.85
        DA2=nan(size(DA2));
        W2=nan(size(W2));
    end
    
    %pairwise information for each combination
    %combinationName, pairwise combo, pairwise combo w/o concentrations, metric, interaction
    N2Info=cell(length(DA2),6);
    ix_2=0;
    for d1=1:nsingles
        for d2=d1+1:nsingles
            ix_2=ix_2+1;
            d2d2=sort([singleNames(d1),singleNames(d2)]); d1d2=strjoin(d2d2,'+');
            N2Info(ix_2,:)={combinationName,d1d2 ,deleteConc(d1d2),DA2(ix_2),W2(ix_2),'not yet determined'};
        end
    end 

%% 3-way interactions
    [~,~,~,~,~,~,~,~,~,E3_123]=calculateE3(w1,w2,w3,w12,w13,w23,w123,CheckInconclusive,rescaleMethod);
    [~,~,~,~,~,~,~,~,~,E3_124]=calculateE3(w1,w2,w4,w12,w14,w24,w124,CheckInconclusive,rescaleMethod);
    [~,~,~,~,~,~,~,~,~,E3_125]=calculateE3(w1,w2,w5,w12,w15,w25,w125,CheckInconclusive,rescaleMethod);
    [~,~,~,~,~,~,~,~,~,E3_134]=calculateE3(w1,w3,w4,w13,w14,w34,w134,CheckInconclusive,rescaleMethod);
    [~,~,~,~,~,~,~,~,~,E3_135]=calculateE3(w1,w3,w5,w13,w15,w35,w135,CheckInconclusive,rescaleMethod);
    [~,~,~,~,~,~,~,~,~,E3_145]=calculateE3(w1,w4,w5,w14,w15,w45,w145,CheckInconclusive,rescaleMethod); %corrected extra XT->XK
    [~,~,~,~,~,~,~,~,~,E3_234]=calculateE3(w2,w3,w4,w23,w24,w34,w234,CheckInconclusive,rescaleMethod);
    [~,~,~,~,~,~,~,~,~,E3_235]=calculateE3(w2,w3,w5,w23,w25,w35,w235,CheckInconclusive,rescaleMethod);
    [~,~,~,~,~,~,~,~,~,E3_245]=calculateE3(w2,w4,w5,w24,w25,w45,w245,CheckInconclusive,rescaleMethod);
    [~,~,~,~,~,~,~,~,~,E3_345]=calculateE3(w3,w4,w5,w34,w35,w45,w345,CheckInconclusive,rescaleMethod);

    [~,DA_123,~,w_123]=calculateN3(w1,w2,w3,w12,w13,w23,w123,CheckInconclusive);
    [~,DA_124,~,w_124]=calculateN3(w1,w2,w4,w12,w14,w24,w124,CheckInconclusive);
    [~,DA_125,~,w_125]=calculateN3(w1,w2,w5,w12,w15,w25,w125,CheckInconclusive);
    [~,DA_134,~,w_134]=calculateN3(w1,w3,w4,w13,w14,w34,w134,CheckInconclusive);
    [~,DA_135,~,w_135]=calculateN3(w1,w3,w5,w13,w15,w35,w135,CheckInconclusive);
    [~,DA_145,~,w_145]=calculateN3(w1,w4,w5,w14,w15,w45,w145,CheckInconclusive);
    [~,DA_234,~,w_234]=calculateN3(w2,w3,w4,w23,w24,w34,w234,CheckInconclusive);
    [~,DA_235,~,w_235]=calculateN3(w2,w3,w5,w23,w25,w35,w235,CheckInconclusive);
    [~,DA_245,~,w_245]=calculateN3(w2,w4,w5,w24,w25,w45,w245,CheckInconclusive);
    [~,DA_345,~,w_345]=calculateN3(w3,w4,w5,w34,w35,w45,w345,CheckInconclusive);
    
    DA3=[DA_123,DA_124,DA_125,DA_134,DA_135,DA_145,DA_234,DA_235,DA_245,DA_345];
    E3=[E3_123,E3_124,E3_125,E3_134,E3_135,E3_145,E3_234,E3_235,E3_245,E3_345];
    W3=[w_123,w_124,w_125,w_134,w_135,w_145,w_234,w_235,w_245,w_345];
    
    if CheckSingleCon && max([w1,w2,w3,w4,w5])>0.85
        DA3=nan(size(DA3)); E3=nan(size(E3)); W3=nan(size(W3));
    end
    
    %triple information for each combination
    %combinationName, triple combo, triple combo w/o concentrations, metric, interaction
    [N3Info,E3Info]=deal(cell(length(DA3),6));
    ix_3=0;
    for d1=1:nsingles
        for d2=d1+1:nsingles
            for d3=d2+1:nsingles
            ix_3=ix_3+1;
            d1d2d3=sort([singleNames(d1),singleNames(d2),singleNames(d3)]); d1d2d3=strjoin(d1d2d3,'+');
            N3Info(ix_3,:)={combinationName,d1d2d3 ,deleteConc(d1d2d3),DA3(ix_3),W3(ix_3),'not yet determined'};
            E3Info(ix_3,:)={combinationName,d1d2d3 ,deleteConc(d1d2d3),E3(ix_3),W3(ix_3),'not yet determined'};
            end
        end
    end

%% 4-way interactions
   
    [DA_1234,E4_1234,w_1234]=calculateN4andE4(w1,w2,w3,w4,w12,w13,w14,w23,w24,w34,w123,w124,w134,w234,w1234,CheckInconclusive,method,rescaleMethod);
    [DA_1235,E4_1235,w_1235]=calculateN4andE4(w1,w2,w3,w5,w12,w13,w15,w23,w25,w35,w123,w125,w135,w235,w1235,CheckInconclusive,method,rescaleMethod);
    [DA_1245,E4_1245,w_1245]=calculateN4andE4(w1,w2,w4,w5,w12,w14,w15,w24,w25,w45,w124,w125,w145,w245,w1245,CheckInconclusive,method,rescaleMethod);
    [DA_1345,E4_1345,w_1345]=calculateN4andE4(w1,w3,w4,w5,w13,w14,w15,w34,w35,w45,w134,w135,w145,w345,w1345,CheckInconclusive,method,rescaleMethod);
    [DA_2345,E4_2345,w_2345]=calculateN4andE4(w2,w3,w4,w5,w23,w24,w25,w34,w35,w45,w234,w235,w245,w345,w2345,CheckInconclusive,method,rescaleMethod);

    %XYZT XYZK XYTK XZTK YZTK
    DA4=[DA_1234,DA_1235,DA_1245,DA_1345,DA_2345];
    E4=[E4_1234,E4_1235,E4_1245,E4_1345,E4_2345];
    W4=[w_1234,w_1235,w_1245,w_1345,w_2345];
    
    if CheckSingleCon && max([w1,w2,w3,w4,w5])>0.85
        DA4=nan(size(DA4)); E4=nan(size(E4));
    end

    %4-drug information for each combination
    %combinationName, 4-drug combo, 4-drug combo w/o concentrations, metric, interaction
    [N4Info,E4Info]=deal(cell(length(DA4),6));
    ix_4=0;
    for d1=1:nsingles
        for d2=d1+1:nsingles
            for d3=d2+1:nsingles
                for d4=d3+1:nsingles
                ix_4=ix_4+1;
                d1d2d3d4=sort([singleNames(d1),singleNames(d2),singleNames(d3),singleNames(d4)]); d1d2d3d4=strjoin(d1d2d3d4,'+');
                N4Info(ix_4,:)={combinationName,d1d2d3d4 ,deleteConc(d1d2d3d4),DA4(ix_4),W4(ix_4),'not yet determined'};
                E4Info(ix_4,:)={combinationName,d1d2d3d4 ,deleteConc(d1d2d3d4),E4(ix_4),W4(ix_4),'not yet determined'};
                end
            end
        end
    end
    
%% 5-way interactions
    %XYZT
    [DA_12345,E5_12345, w_12345]=calculateN5andE5(w1,w2,w3,w4,w5,...
                                w12,w13,w14,w15,w23,w24,w25,w34,w35,w45,...
                                w123,w124,w125,w134,w135,w145,w234,w235,w245,w345,...
                                w1234,w1235,w1245,w1345,w2345,...
                                w12345,...
                                CheckInconclusive,method,rescaleMethod);
                            
    %XYZTK 
    DA5=[DA_12345];
    E5=[E5_12345];
    W5=[w_12345];
    
    if CheckSingleCon && max([w1,w2,w3,w4,w5])>0.85
        DA5=nan(size(DA5)); E5=nan(size(E5)); W5=nan(size(W5));  
    end
           
    %5-drug information for each combination
    %combinationName, 5-drug combo, 5-drug combo w/o concentrations, metric, interaction
    [N5Info,E5Info]=deal(cell(length(DA5),6));
    ix_5=0;
    for d1=1:nsingles
        for d2=d1+1:nsingles
            for d3=d2+1:nsingles
                for d4=d3+1:nsingles
                    for d5=d4+1:nsingles
                    ix_5=ix_5+1;
                    d1d2d3d4d5=sort([singleNames(d1),singleNames(d2),singleNames(d3),singleNames(d4),singleNames(d5)]); d1d2d3d4d5=strjoin(d1d2d3d4d5,'+');
                    N5Info(ix_5,:)={combinationName,d1d2d3d4d5 ,deleteConc(d1d2d3d4d5),DA5(ix_5),W5(ix_5),'not yet determined'};
                    E5Info(ix_5,:)={combinationName,d1d2d3d4d5 ,deleteConc(d1d2d3d4d5),E5(ix_5),W5(ix_5),'not yet determined'};
                    end
                end
            end
        end
    end
    
    metrics=num2cell([DA2,DA3,E3,DA4,E4,DA5,E5]);
end