%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function circint(xc1,yc1,r1,xc2,yc2,r2,numbers,percentage_numbers)
%     numbers=[1 1 1]
%     xc1=0; yc1=0; r1=1;
%     xc2=1.1; yc2=0; r2=1;
    % xc1, yc1, r1: x,y coordinates of center and radius of first circle
    % xc2, yc2, r2: x,y coordinates of center and radius of second circle

    if isnumeric(numbers), numbers=num2cell(numbers); end 
    
    if nargin<=7
        strDA = {'Net','','Interaction','','',numbers{1}};
        strboth = {'Both','','',numbers{2}};
        strE3 = {'Emergent','','Interaction','','',numbers{3}};
    else 
        if isnumeric(percentage_numbers), percentage_numbers=num2cell(percentage_numbers); end 

        strDA = {'Net','Interaction','','',numbers{1},percentage_numbers{1}};
        strboth = {'Both','',numbers{2},percentage_numbers{2}};
        strE3 = {'Emergent','Interaction','','',numbers{3},percentage_numbers{3}};
    end
    %Polar coordinates
    theta = 0:2*pi/100:2*pi;

    %Circle equations
    x1 = r1*cos(theta) + xc1; y1 = r1*sin(theta) + yc1;

    x2 = r2*cos(theta) + xc2; y2 = r2*sin(theta) + yc2;

    %find minimum & maximum values of x, y
    xmin = min([xc1-r1,xc2-r2]); xmax = max([xc1+r1,xc2+r2]);

    ymin = min([yc1-r1,yc2-r2]); ymax = max([yc1+r1,yc2+r2]);

    %Number of grid cells (increase this value to get more accuracy in shading)
    ndiv = 1000;

    %Create mesh
    xa = xmin:(xmax-xmin)/ndiv:xmax;
    ya = ymin:(ymax-ymin)/ndiv:ymax; 
    [xo,yo] = meshgrid(xa,ya);

    %Find Intersection Region
    G = double( ((xo - xc1).^2 + (yo - yc1).^2 < r1^2) & ((xo - xc2).^2 + (yo -yc2).^2 < r2^2));

    %Replace zeros with NaN (Better look with pcolor)
    G(G == 0) = NaN;

    %Show plots
    
%     N=1000;
%     filledCircle([xc1,yc1],r1,N,'r')
    
    h=fill(xa,ya,'r');
    pcolor(xa,ya,G)
    colormap(gray)
    shading interp
    grid off
    hold on
    
     plot(x1,y1,'k')
    axis equal
    hold on
    plot(x2,y2,'k')
    axis tight
    
    

dim = [-.5 .2 .3 .3];
    
% % annotation('textbox',dim,'String',str,'FitBoxToText','on');
% text(-(r1-xc1)/2,(r1-yc1)/2,str,'HorizontalAlignment','Center','FontName','Arial','FontSize',16)

text(-0.3,0.2,strDA,'HorizontalAlignment','Center','fontsize',24,'FontName','Arial')
text(0.5,0.2,strboth,'HorizontalAlignment','Center','fontsize',24,'FontName','Arial')
text(1.4,0.2,strE3,'HorizontalAlignment','Center','fontsize',24,'FontName','Arial')

ylim([ymin,ymax+0.2]);

set(gca,'fontsize',24,'FontName','Arial'); axis off;
end 
    