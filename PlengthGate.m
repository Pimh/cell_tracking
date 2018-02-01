
function [valid_xPlength,valid_yPlength,valid_Plength,finalX,finalY]=PlengthGate(fname,refFile,type,timeStep)
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
%xdisp=[]; ydisp=[];
valid_xvel=[]; valid_yvel=[];
valid_stdx=[]; valid_stdy=[];
pixpermic= 1/1;
Spd_cutoff = 0.2;
P_cutoff=0.2;
%x_cutoff = 2;
%y_cutoff =3;
i=0; mint=inf;

fid = fopen(['GatedPlength' fname '.txt'],'w');
fprintf(fid,'\n%s',fname);
fprintf(fid,'speed threshold:      %d\n', Spd_cutoff);
fprintf(fid,'  Speed        x-velocity       y-velocity       x-Plength       y-Plength      Plength\n\n');

colrVec = colorInterp(nc);
for ind=1:nc
    [x,y,t]=getCell(cells,ind);
    x = x.*pixpermic;
    y = y.*pixpermic;
    t = t.*timeStep;
    
    if t(length(t))< mint
        mint = t(length(t));
    end
    
    allSpd(ind) = calSpeed(x,y,t);
    %[xvel(ind),yvel(ind)] = calVelocity(x,y,t);
    %[xPlength(ind),yPlength(ind)]= calPlength(x,y);
    [stdx(ind),stdy(ind)]= calStd(x,y);
    [xPlength(ind),yPlength(ind),Plength(ind)]= calPlength(x,y,t);
    [xvel(ind),yvel(ind)] = calVelocity(x,y,t);
    
    if  Plength(ind)>P_cutoff
        i=i+1;
        valid_Spd = [valid_Spd allSpd(ind)];
        valid_xvel= [valid_xvel xvel(ind)];
        valid_yvel= [valid_yvel yvel(ind)];
        valid_xPlength= [valid_xPlength xPlength(ind)];
        valid_yPlength= [valid_yPlength yPlength(ind)];
        valid_Plength= [valid_Plength Plength(ind)];
        %valid_motility= [valid_motility motility(ind)];
        
        finalX(i)=x(length(x))-x(1);
        finalY(i)=y(length(y))-y(1);
      
        fprintf(fid,'%10.6f', valid_Spd(i));
%         fprintf(fid,'     %10.6f', valid_stdx(i))
%         fprintf(fid,'     %10.6f', valid_stdy(i))
        fprintf(fid,'     %10.6f', valid_xvel(i));
        fprintf(fid,'     %10.6f', valid_yvel(i));
        fprintf(fid,'     %10.6f', valid_xPlength(i));
        fprintf(fid,'     %10.6f', valid_yPlength(i));
        fprintf(fid,'     %10.6f\n', valid_Plength(i));
        %fprintf(fid,'     %10.6f\n', valid_motility(i));
        
        
    end
    %Figure1: Trajectory plot
    hold on
    axis equal
    axis([-200 200 -200 200])
    grid on
    figure(1)
    colr = colrVec(ind,:);
    plot(x-x(1),y-y(1),'LineStyle','-','Color',colr)
    cell=sprintf('%d',ind);
    %title(cell)
    %pause
end


%Figure2: Histogram of speed
figure(2)
hist(valid_Spd,10)
title('speed');
avgSpd = mean(valid_Spd);
stdSpd = std(valid_Spd);

%Figure3: Histogram of x-velocity
figure(3)
hist(valid_xvel,10)
title('x-velocity');
avgXVel = mean(valid_xvel);
stdXVel = std(valid_xvel);

%Figure4: Histogram of y-velocity
figure(4)
hist(valid_yvel,10)
title('y-velocity');
avgYVel = mean(valid_yvel);
stdYVel = std(valid_yvel);


% %Figure5: Histogram of stdX
% figure(5)
% hist(valid_stdx,20)
% title('stdX');
% 
% %Figure6: Histogram of stdY
% figure(6)
% hist(valid_stdy,20)
% title('stdY');

avgXPlength = mean(valid_xPlength);
stdXPlength = std(valid_xPlength);
avgYPlength = mean(valid_yPlength);
stdYPlength = std(valid_yPlength);
avgPlength = mean(valid_Plength);
stdPlength = std(valid_Plength);
% avgMotility = mean(valid_motility);
% stdMotility = std(valid_motility);

% avgstdx = mean(valid_stdx);
% stdStdx = std(valid_stdx);
% avgstdy = mean(valid_stdy);
% stdStdy = std(valid_stdy);
%fid = fopen('Results.txt','a');

fprintf(fid,'\n%s\n',fname);
fprintf(fid,'threshold:      %d\n', Spd_cutoff);
fprintf(fid,'# total cells:      %d\n', nc);
fprintf(fid,'# motile cells:     %d\n', i);
fprintf(fid,'--------------------- Average --------  Stadev ---------------------\n\n');
fprintf(fid,'speed:              %10.7f     %10.10f\n', avgSpd,stdSpd);
fprintf(fid,'x-velocity:         %10.7f     %10.10f\n', avgXVel,stdXVel);
fprintf(fid,'y-velocity:         %10.7f     %10.10f\n', avgYVel,stdYVel);
fprintf(fid,'x-persistentlength: %10.7f     %10.10f\n', avgXPlength,stdXPlength);
fprintf(fid,'y-persistentlength: %10.7f     %10.10f\n', avgYPlength,stdYPlength);
fprintf(fid,'persistentlength  : %10.7f     %10.10f\n', avgPlength,stdPlength);
%fprintf(fid,'motility  :         %10.7f     %10.10f\n', avgMotility,stdMotility);
% fprintf(fid,'stdx:                    %10.10f\n', mean(valid_stdx))
% fprintf(fid,'stdy:                    %10.10f\n', mean(valid_stdy))
fprintf(fid,'mininum time duration:   %10d\n', mint);

fclose(fid);

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
function [xPlength,yPlength,Plength]= calPlength(x,y,t)

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
Fcell=fopen('cell.txt','w');
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
fclose(Fcell);   