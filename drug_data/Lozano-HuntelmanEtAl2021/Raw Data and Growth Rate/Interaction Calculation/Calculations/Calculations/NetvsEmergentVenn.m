function NetvsEmergentVenn(vennstatistics,titletext,vennnumbers)

%vennstatistics: [onlyDA, both, onlyE3]
    xc1=0; yc1=0; r1=1;
    xc2=1.1; yc2=0; r2=1;
    if nargin<=2
        circint(xc1,yc1,r1,xc2,yc2,r2,vennstatistics)
    else 
        circint(xc1,yc1,r1,xc2,yc2,r2,vennstatistics,vennnumbers)
    end
    title(titletext,'fontsize',24,'FontName','Arial','FontWeight','normal');

end