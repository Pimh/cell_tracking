function xdis= poissonDist()
%Read excel file
%dataset2,4,6


[cells1,nc(1)]= CalibratedPos('D3P1.xls','D3P1ref.xls',3);
[cells2,nc(2)]= CalibratedPos('D3P2.xls','D3P2ref.xls',3);

pixpermic=1;
xdis= cell(120,1);

for set=1:2
    
    if set ==1
        cells=cells1;
    else
        cells=cells2;
    end
    for ind=1:nc
        x=[];y=[];t=[];deltax=[];
        [x,y,t]=getCell(cells,ind);
        x = x.*pixpermic;
        y = y.*pixpermic;
        t = t-t(1)+1; 
        deltax= x-x(1);


        for i=2:length(t)
            xdis{t(i)}(length(xdis{t(i)})+1) =  deltax(i);               
        end   

    end
end

n24=hist(xdis{24},x);
n48=hist(xdis{48},x);
n72=hist(xdis{72},x);
n96=hist(xdis{96},x);
n120=hist(xdis{120},x);
plot(x,n24/length(xdis{24}),x,n48/length(xdis{48}),x,n72/length(xdis{72}),...
    x,n96/length(xdis{96}),x,n120/length(xdis{120}))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% sorting out x-coordinates, y-coordinates, and time vector %%%%%%%%%%%
function [x,y,t]=getCell(cells,target)

x = cells{target}(:,1);
y = cells{target}(:,2);
t = cells{target}(:,3);

