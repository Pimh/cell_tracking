
function phenotype(fname,refFile,type,timeStep)
%Include reference file
%Display info of every single cell that is qualified as a motile cell
%Result is printed on the different files for different device positions

if type==3
    xc=1;yc=2;tc=4;IDc=5;
elseif type ==2
    xc=4;yc=5;tc=3;IDc=2;
end
[cells,nc]=loadCellData([fname '.txt'],[refFile '.txt'],xc,yc,tc,IDc);

cellSpd=[]; valid_Spd=[]; VUratio=[];
xvel=[]; yvel=[];
xPlength=[]; yPlength=[];
stdx=[]; stdy=[];
%xdisp=[]; ydisp=[];
valid_xvel=[]; valid_yvel=[];
valid_stdx=[]; valid_stdy=[];
pixpermic= 1/1;
%Spd_cutoff = 0.1;
x_cutoff = 3;
y_cutoff =3;
i=0; mint=inf;


colrVec = colorInterp(nc);

a=0;b=0;Spd_A=[];Spd_B=[];xPl_A=[];xPl_B=[];Pl_A=[];Pl_B=[];Vx_A=[];Vx_B=[];
for ind=1:nc
    [x,y,t]=getCell(cells,ind);
    x = x.*pixpermic;
    y = y.*pixpermic;
    t = t.*timeStep;
    
    if t(length(t))< mint
        mint = t(length(t));
    end
    
    [allSpd(ind),xSpeed(ind)] = calSpeed(x,y,t);
    %[xvel(ind),yvel(ind)] = calVelocity(x,y,t);
    %[xPlength(ind),yPlength(ind)]= calPlength(x,y);
    [stdx(ind),stdy(ind)]= calStd(x,y);
    [xPlength(ind),yPlength(ind),Plength(ind)]= calPlength(x,y);
    [xvel(ind),yvel(ind)] = calVelocity(x,y,t);
    VUratio = [VUratio xvel(ind)/allSpd(ind)];
    if abs(xvel(ind)/allSpd(ind))>0.05
        a=a+1;
        Spd_A= [Spd_A allSpd(ind)];
        Vx_A = [Vx_A xvel(ind)];
        xPl_A= [xPl_A xPlength(ind)];
        Pl_A = [Pl_A sqrt(xPlength(ind)^2+yPlength(ind)^2)];
    else
        b=b+1;
        Spd_B= [Spd_B allSpd(ind)];
        Vx_B = [Vx_B xvel(ind)];
        xPl_B= [xPl_B xPlength(ind)];
        Pl_B = [Pl_B sqrt(xPlength(ind)^2+yPlength(ind)^2)];
    end
    
%     hold on
%     figure(1)
%     plot(t/timeStep,x-x(1),'color',colrVec(ind,:))
%     grid on
%     axis([0 200 -100 100])
%     xlabel('time point')
%     ylabel('x displacement')
%     title(fname)
%     
%     hold on 
%     figure(2)
%     plot(t/timeStep,y-y(1),'color',colrVec(ind,:))
%     grid on
%     axis([0 200 -100 100])
%     xlabel('time point')
%     ylabel('y displacement')
%     title(fname)
    
%     hold on
%     figure(3)
%     plot(allSpd(ind),xvel(ind),'k*')
%     %axis([0 0.6 -0.6 .6])
%     xlabel('speed')
%     ylabel('x velocity')
% 
%     hold on
%     figure(4)
%     plot(xSpeed(ind),xvel(ind),'k*')
%     axis([0 0.35 -.35 .35])
%     xlabel('x speed')
%     ylabel('x velocity')

%     hold on
%     figure(5)
%     plot(allSpd(ind),yvel(ind),'k*')
%     %axis([0 0.6 -0.6 .6])
%     xlabel('speed')
%     ylabel('y velocity')

      hold on
    figure(6)
    plot(Plength(ind),xPlength(ind),'k*')
    axis equal
    xlabel('x Persistent length')
    ylabel('Persistent length')
    
    %pause
end
a
AvgU_A=mean(Spd_A)
AvgVx_A=mean(Vx_A)
AvgPx_A=mean(xPl_A)
AvgP_A=mean(Pl_A)

b
AvgU_B=mean(Spd_B)
AvgVx_B=mean(Vx_B)
AvgPx_B=mean(xPl_B)
AvgP_B=mean(Pl_B)

 figure(1)
% hist(VUratio)
%sum(VUratio>0.05)
hist(xPlength./Plength)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% calculating for speed %%%%%%%%%%%%%%%%%%%%%%%%
function [cellSpd,xCellSpd] = calSpeed(x,y,t)

L= length(x);cellSpd=0;xSpd=0;
for i=2:L
    xdist = x(i)-x(i-1);
    ydist = y(i)-y(i-1);
    sdist = sqrt(xdist^2+ydist^2);
    tInt = t(i)-t(i-1);
    spd = sdist/(tInt);
    
    xSpd = xSpd + abs(xdist)/tInt;
    cellSpd = cellSpd + spd;
end

xCellSpd= xSpd/L;
cellSpd = cellSpd/L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% calculating for velocity %%%%%%%%%%%%%%%%%%%%%%%%
function [xvel,yvel] = calVelocity(x,y,t)

L= length(x);sumx=0;sumy=0;
for i=2:L
    xdist = x(i)-x(i-1);
    sumx = sumx + xdist;
    ydist = y(i)-y(i-1);
    sumy = sumy + ydist;
end
xvel = sumx/(t(L));
yvel = sumy/(t(L));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% calculating for persistent length %%%%%%%%%%%%%%%%%%%%%%%%
function [xPlength,yPlength,Plength]= calPlength(x,y)

L=length(x);sumdist=0;
for i=2:L
    xdist = x(i)-x(i-1);
    ydist = y(i)-y(i-1);
    sdist = sqrt(xdist^2+ydist^2);
    
    sumdist = sumdist+sdist;
end
xdisp= x(L)-x(1);
ydisp= y(L)-y(1);
Plength= sqrt(xdisp^2+ydisp^2)/sumdist;
xPlength= xdisp/sumdist;
yPlength= ydisp/sumdist;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% calculating for persistent length
function [stdx,stdy]= calStd(x,y)

stdx= std(x);
stdy= std(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% sorting out x-coordinates, y-coordinates, and time vector %%%%%%%%%%%
function [x,y,t]=getCell(cells,target)

x = cells{target}(:,1);
y = cells{target}(:,2);
t = cells{target}(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% read the raw data and store the read value into separate element of
%%%%%%%%%%%% cell array cells
function [cells,nc]=loadCellData(fname,refFile,xc,yc,tc,IDc)

particles= load(fname);
[r,c] = size(particles);

particle= [particles(:,xc) particles(:,yc) particles(:,tc) particles(:,IDc)];

nc=0; cells={};i=1; 
% cont=1;
while i <=r 
%     if particle(i,3) == 1
      %does not have to start at tp=1  
        first =i;
        i=i+1;
        while i <= r && (i ==1 ||particle(i,4)==particle(i-1,4))
            i=i+1;
        end
        
        if i ==r
            last = r;
        else
            last = i-1;
        end
        nc = nc+1;
        cells{nc} = particle(first:last,:);
%     else 
%         cont=0;
%     end
end

%%%%%%%reference file%%%%%%%
Ref = load(refFile);

[n,m] = size(Ref);
xRef=[]; yRef=[];

for k= 1:n
    xRef(k)= Ref(k,1)-Ref(1,1);
    yRef(k)= Ref(k,2)-Ref(1,2);
end
Fcell=fopen('cell.txt','a')
for j = 1:nc
    fprintf(Fcell,'\n\n%d\n',j);
    L= length(cells{j}(:,1));
    
%     if length(xRef) == L
%         cells{j}(:,1) = cells{j}(:,1)-(xRef(1:L))';
%         cells{j}(:,2) = cells{j}(:,2)-(yRef(1:L))';
%         
%     else
        c=1;r=1;
        while r<=length(xRef) && c<=L
            
            if c>1 && cells{j}(c-1,3)>cells{j}(c,3)
                fprintf('\ncell# = %d', j);
            end
            
            if Ref(r,tc)== cells{j}(c,3)
                cells{j}(c,1) = cells{j}(c,1)-xRef(r);
                cells{j}(c,2) = cells{j}(c,2)-yRef(r);
                
                fprintf(Fcell,'%10d',cells{j}(c,3));
                fprintf(Fcell,'          %d',cells{j}(c,4));
                fprintf(Fcell,'          %10.2f',cells{j}(c,1));
                fprintf(Fcell,'          %10.2f\n',cells{j}(c,2));
                c=c+1;
                r=r+1;
            else
                r=r+1;
            end
            
        end
%     end
    
    
    
end
fclose(Fcell)   