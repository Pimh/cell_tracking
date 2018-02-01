function [valid_Spd1,valid_xvel1,valid_Plength1]=PXversusP(fname,refFile,type,timeStep)
%Include reference file
%Display info of every single cell that is qualified as a motile cell
%Result is printed on the different files for different device positions

if type==3
    xc=1;yc=2;tc=4;IDc=5;
elseif type ==2
    xc=4;yc=5;tc=3;IDc=2;
end
[cells,nc]=loadCellData([fname '.txt'],[refFile '.txt'],xc,yc,tc,IDc);

cellSpd=[]; valid_Spd=[];
xvel=[]; yvel=[];
valid_xPlength=[]; valid_yPlength=[];
valid_Plength=[]; valid_motility=[];
stdx=[]; stdy=[];
valid_xvel=[]; valid_yvel=[];
valid_stdx=[]; valid_stdy=[];
pixpermic= 1/1;
threshold=0.2;
i=0; mint=inf;j=0;

fid1 = fopen(['A_' fname '.txt'],'w');
fprintf(fid1,'\n%s',fname)
fprintf(fid1,'\n %2.1f  \n\n',threshold)
fprintf(fid1,'#  Speed        x-velocity       y-velocity       x-Plength       y-Plength      Plength       motility\n\n')

fid2 = fopen(['B_' fname '.txt'],'w');
fprintf(fid2,'\n%s',fname)
fprintf(fid2,'\n %2.1f  \n\n',threshold)
fprintf(fid2,'#  Speed        x-velocity       y-velocity       x-Plength       y-Plength      Plength       motility\n\n')

colrVec = colorInterp(nc);
for ind=1:nc
    [x,y,t]=getCell(cells,ind);
    x = x.*pixpermic;
    y = y.*pixpermic;
    t = t.*timeStep;
    
%     if t(length(t))< mint
%         mint = t(length(t));
%     end
    
    allSpd(ind) = calSpeed(x,y,t);
    [stdx(ind),stdy(ind)]= calStd(x,y);
    [xPlength(ind),yPlength(ind),Plength(ind),motility(ind)]= calPlength(x,y,t);
    [xvel(ind),yvel(ind)] = calVelocity(x,y,t);
    
    if Plength(ind)>threshold
        i=i+1;
        valid_Spd1(i) =  allSpd(ind);
        valid_xvel1(i)=  xvel(ind);
        valid_yvel1(i)=  yvel(ind);
        valid_xPlength1(i)=  xPlength(ind);
        valid_yPlength1(i)=  yPlength(ind);
        valid_Plength1(i)=  Plength(ind);
        valid_motility1(i)=  motility(ind);
        
        
      
        fprintf(fid1,'%3d', ind)
        fprintf(fid1,'     %10.6f', valid_Spd1(i))
        fprintf(fid1,'     %10.6f', valid_xvel1(i))
        fprintf(fid1,'     %10.6f', valid_yvel1(i))
        fprintf(fid1,'     %10.6f', valid_xPlength1(i))
        fprintf(fid1,'     %10.6f', valid_yPlength1(i))
        fprintf(fid1,'     %10.6f', valid_Plength1(i))
        fprintf(fid1,'     %10.6f\n', valid_motility1(i))
    else
        j=j+1;
        valid_Spd2(j) =  allSpd(ind);
        valid_xvel2(j)=  xvel(ind);
        valid_yvel2(j)=  yvel(ind);
        valid_xPlength2(j)=  xPlength(ind);
        valid_yPlength2(j)=  yPlength(ind);
        valid_Plength2(j)=  Plength(ind);
        valid_motility2(j)=  motility(ind);
        
        
      
        fprintf(fid2,'%3d', ind)
        fprintf(fid2,'     %10.6f', valid_Spd2(j))
        fprintf(fid2,'     %10.6f', valid_xvel2(j))
        fprintf(fid2,'     %10.6f', valid_yvel2(j))
        fprintf(fid2,'     %10.6f', valid_xPlength2(j))
        fprintf(fid2,'     %10.6f', valid_yPlength2(j))
        fprintf(fid2,'     %10.6f', valid_Plength2(j))
        fprintf(fid2,'     %10.6f\n', valid_motility2(j))
        
        
        
    end
    %Figure1: Trajectory plot
%     hold on
%     axis equal
%     axis([-200 200 -200 200])
%     grid on
%     figure(1)
%     colr = colrVec(ind,:);
%     plot(x-x(1),y-y(1),'LineStyle','-','Color',colr)
%     cell=sprintf('%d',ind);
%     %title(cell)
%     %pause
    
    hold on
    figure(1)
    plot(Plength(ind),xPlength(ind),'ro','MarkerFace','r','MarkerSize',5)
    axis([0 1 -1 1])
    xlabel('Persistence length','FontSize',14)
    ylabel('x Persistence length','FontSize',14)
    
    
    
    
end
% hold on
% figure(1)
% hist(xPlength./Plength)
% title('distribution of Px/P')
% no=sum(abs(xPlength./Plength)>0.6)

avgSpd1 = mean(valid_Spd1);
stdSpd1 = std(valid_Spd1);

avgXVel1 = mean(valid_xvel1);
stdXVel1 = std(valid_xvel1);

avgYVel1 = mean(valid_yvel1);
stdYVel1 = std(valid_yvel1);

avgXPlength1 = mean(valid_xPlength1);
stdXPlength1 = std(valid_xPlength1);

avgYPlength1 = mean(valid_yPlength1);
stdYPlength1 = std(valid_yPlength1);

avgPlength1 = mean(valid_Plength1);
stdPlength1 = std(valid_Plength1);

avgMotility1 = mean(valid_motility1);
stdMotility1 = std(valid_motility1);


fprintf(fid1,'\n%s\n',fname)
fprintf(fid1,'# total cells:      %d\n', nc)
fprintf(fid1,'# motile cells:     %d\n', i)
fprintf(fid1,'--------------------- Average --------  Stadev ---------------------\n\n')
fprintf(fid1,'speed:              %10.7f     %10.10f\n', avgSpd1,stdSpd1)
fprintf(fid1,'x-velocity:         %10.7f     %10.10f\n', avgXVel1,stdXVel1)
fprintf(fid1,'y-velocity:         %10.7f     %10.10f\n', avgYVel1,stdYVel1)
fprintf(fid1,'x-persistentlength: %10.7f     %10.10f\n', avgXPlength1,stdXPlength1)
fprintf(fid1,'y-persistentlength: %10.7f     %10.10f\n', avgYPlength1,stdYPlength1)
fprintf(fid1,'persistentlength  : %10.7f     %10.10f\n', avgPlength1,stdPlength1)
fprintf(fid1,'motility  :         %10.7f     %10.10f\n', avgMotility1,stdMotility1)

avgSpd2 = mean(valid_Spd2);
stdSpd2 = std(valid_Spd2);

avgXVel2 = mean(valid_xvel2);
stdXVel2 = std(valid_xvel2);

avgYVel2 = mean(valid_yvel2);
stdYVel2 = std(valid_yvel2);

avgXPlength2 = mean(valid_xPlength2);
stdXPlength2 = std(valid_xPlength2);

avgYPlength2 = mean(valid_yPlength2);
stdYPlength2 = std(valid_yPlength2);

avgPlength2 = mean(valid_Plength2);
stdPlength2 = std(valid_Plength2);

avgMotility2 = mean(valid_motility2);
stdMotility2 = std(valid_motility2);


fprintf(fid2,'\n%s\n',fname)
fprintf(fid2,'# total cells:      %d\n', nc)
fprintf(fid2,'# motile cells:     %d\n', j)
fprintf(fid2,'--------------------- Average --------  Stadev ---------------------\n\n')
fprintf(fid2,'speed:              %10.7f     %10.10f\n', avgSpd2,stdSpd2)
fprintf(fid2,'x-velocity:         %10.7f     %10.10f\n', avgXVel2,stdXVel2)
fprintf(fid2,'y-velocity:         %10.7f     %10.10f\n', avgYVel2,stdYVel2)
fprintf(fid2,'x-persistentlength: %10.7f     %10.10f\n', avgXPlength2,stdXPlength2)
fprintf(fid2,'y-persistentlength: %10.7f     %10.10f\n', avgYPlength2,stdYPlength2)
fprintf(fid2,'persistentlength  : %10.7f     %10.10f\n', avgPlength2,stdPlength2)
fprintf(fid2,'motility  :         %10.7f     %10.10f\n', avgMotility2,stdMotility2)

fclose(fid1);
fclose(fid2);

%xlswrite('test.xls',valid_Spd','speed','A1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% calculating for speed %%%%%%%%%%%%%%%%%%%%%%%%
function cellSpd = calSpeed(x,y,t)

L= length(x);cellSpd=0;
for i=2:L
    xdist = x(i)-x(i-1);
    ydist = y(i)-y(i-1);
    sdist = sqrt(xdist^2+ydist^2);
    tInt = t(i)-t(i-1);
    spd = sdist/(tInt);
    
    cellSpd = cellSpd + spd;
end

cellSpd = cellSpd/(L-1);

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
xvel = sumx/(t(L)-t(1));
yvel = sumy/(t(L)-t(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% calculating for persistent length %%%%%%%%%%%%%%%%%%%%%%%%
function [xPlength,yPlength,Plength,motility]= calPlength(x,y,t)

L=length(x);sumdist=0;sumdist_square=0;
for i=2:L
    xdist = x(i)-x(i-1);
    ydist = y(i)-y(i-1);
    dist_square= xdist^2+ydist^2;
    sdist = sqrt(dist_square);
    
    sumdist = sumdist+sdist;
    sumdist_square= sumdist_square+dist_square;
end
xdisp= x(L)-x(1);
ydisp= y(L)-y(1);

xPlength= xdisp/sumdist;
yPlength= ydisp/sumdist;
Plength= sqrt(xdisp^2+ydisp^2)/sumdist;
motility= sumdist_square/(t(L)-t(1));

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
Fcell=fopen('cell.txt','w')
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