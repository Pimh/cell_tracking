
function [U,Vx,Vy,P,Px]=stdGate(fname,refFile,type,timeStep)
%Calculate and return migration parameters U,Vx,Vy,P,Px
%Generate a polar plot and histogram of U,Vx,Vy
%Creat an output file,a text file with a name of fname
%Require 4 input parameters: fname,refFile,type, and timeStep
%fname is a string vector of the file name that will be analyzed
%refFile is the reference file used for correct frame drifting
%type = 2 if used 2D tracking and 3 if used 3D/Imaris tracking
%timestep = time interval between 2 time frames specified in minute

if type==3
    xc=1;yc=2;tc=7;IDc=8;
elseif type ==2
    xc=4;yc=5;tc=3;IDc=2;
end

[cells,nc]=loadCellData(fname,refFile,xc,yc,tc,IDc);

%Initialize variables & vectors
U=[];xvel=[]; yvel=[];
Px=[]; Py=[];P=[]; 
stdx=[]; stdy=[];
Vx=[]; Vy=[];
i=0; mint=inf;

%Conversion factor, convert pixel to micron
pixpermic= 1;

%Criteria for motile cells
% stdx and stdy are used as criteria in this case
% any cells with stdx and stdy smaller than the x_cutoff and y_cutoff
% will be not considered in the analysis, but still on the polar plots
x_cutoff = 2;
y_cutoff =2.5;


% Creat a text file for recording migration parameters of an individaul cell
% and also the average values over the cell population of the same data set
fid = fopen(['result' fname '.txt'],'w');
fprintf(fid,'\n%s',fname);
fprintf(fid,'\n x= %2.1f  y= %2.1f\n\n',x_cutoff,y_cutoff);
fprintf(fid,'ID  Speed        x-velocity       y-velocity       x-Plength       y-Plength      Plength    Time Duration\n\n');

%Color vector for polar plot
colrVec = colorInterp(nc);

%Loop over entire data set to calculate for migration parameters of each
%cell
for ind=1:nc
    [x,y,t]=getCell(cells,ind);
    x = x.*pixpermic;
    y = y.*pixpermic;
    t = t.*timeStep;
    
    %Calculate speed, std, persistence length, and velocity
    allSpd(ind) = calSpeed(x,y,t);
    [stdx(ind),stdy(ind)]= calStd(x,y);
    [xPlength(ind),yPlength(ind),Plength(ind)]= calPlength(x,y,t);
    [xvel(ind),yvel(ind)] = calVelocity(x,y,t);
    
    %Screen out non-motile cells
    if stdx(ind) >= x_cutoff && stdy(ind)>= y_cutoff
        i=i+1;
        U = [U allSpd(ind)];
        Vx= [Vx xvel(ind)];
        Vy= [Vy yvel(ind)];
        Px= [Px xPlength(ind)];
        Py= [Py yPlength(ind)];
        P= [P Plength(ind)];
      
        tl= t(length(t))-t(1);
        if tl < mint
            mint = tl;
        end
        
        fprintf(fid,'%d', ind);
        fprintf(fid,'     %10.6f', U(i));
        fprintf(fid,'     %10.6f', Vx(i));
        fprintf(fid,'     %10.6f', Vy(i));
        fprintf(fid,'     %10.6f', Px(i));
        fprintf(fid,'     %10.6f', Py(i));
        fprintf(fid,'     %10.6f', P(i));
        fprintf(fid,'     %d\n', tl);
       
        
        
    end
    
    %Figure1: Polar Plot
    hold on
    axis equal
    axis([-200 200 -200 200])
    grid on
    figure(1)
    colr = colrVec(ind,:);
    plot(x-x(1),y-y(1),'LineStyle','-','LineWidth',2,'Color',colr)
    hold on
    plot(x(length(x))-x(1),y(length(y))-y(1),'ko','MarkerEdgeColor','k',...
                'MarkerFaceColor',[0 0 0],...
                'MarkerSize',5) 
   
end


%Figure2: Histogram of speed
figure(2)
hist(U,10)
title('speed');
avgSpd = mean(U);
stdSpd = std(U);

%Figure3: Histogram of x-velocity
figure(3)
hist(Vx,10)
title('x-velocity');
avgXVel = mean(Vx);
stdXVel = std(Vx);

%Figure4: Histogram of y-velocity
figure(4)
hist(Vy,10)
title('y-velocity');
avgYVel = mean(Vy);
stdYVel = std(Vy);

avgXPlength = mean(Px);
stdXPlength = std(Px);
avgYPlength = mean(Py);
stdYPlength = std(Py);
avgPlength = mean(P);
stdPlength = std(P);


fprintf(fid,'\n%s\n',fname);
fprintf(fid,'# total cells:      %d\n', nc);
fprintf(fid,'# motile cells:     %d\n', i);
fprintf(fid,'--------------------- Average --------  Stadev ---------------------\n\n');
fprintf(fid,'speed:              %10.7f     %10.10f\n', avgSpd,stdSpd);
fprintf(fid,'x-velocity:         %10.7f     %10.10f\n', avgXVel,stdXVel);
fprintf(fid,'y-velocity:         %10.7f     %10.10f\n', avgYVel,stdYVel);
fprintf(fid,'x-persistentlength: %10.7f     %10.10f\n', avgXPlength,stdXPlength);
fprintf(fid,'y-persistentlength: %10.7f     %10.10f\n', avgYPlength,stdYPlength);
fprintf(fid,'persistentlength  : %10.7f     %10.10f\n', avgPlength,stdPlength);
fprintf(fid,'mininum time duration:   %10d\n', mint);

fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% calculate speed %%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%% calculate velocity %%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%% calculate persistence length %%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%% calculate persistence length
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

particles= xlsread(fname);
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
Ref = xlsread(refFile);

[n,m] = size(Ref);
xRef=[]; yRef=[];

for k= 1:n
    xRef(k)= Ref(k,xc)-Ref(1,xc);
    yRef(k)= Ref(k,yc)-Ref(1,yc);
end
Fcell=fopen('cell.txt','w');
for j = 1:nc
    fprintf(Fcell,'\n\n%d\n',j);
    L= length(cells{j}(:,1));

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

end
fclose(Fcell);   