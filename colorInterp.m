
function colrVec = colorInterp(TTtp)
%color interpolation
%TTtp = 20;
blue = [0 0 1];
purple=[1 0 1];
cyan  =[0 1 1];
red   =[1 0 0];
green =[0 1 0];
yellow = [1 1 0];

colrCell= {red yellow green cyan purple blue};
L = length(colrCell)-1;
int= ceil(TTtp/L);
tp=1;
colrVec=[];

for i=1:L
    for j= 1:int
        colr1 = (((j-1).*colrCell{i+1})+ ((int-j).*colrCell{i}))/(int-1);
        colrVec = [colrVec; colr1];
        %plot(x(tp:tp+1),y(tp:tp+1),'LineStyle','-','Color',colr1)
        %x=[1 2 2 1]*tp;
        %fill(x,[1 1 2 2],colr1)
        %pause
        tp = tp+1;
    end
end
