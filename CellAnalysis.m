function [U,Vx,Vy,P,Px,Py]=CellAnalysis(data,fname)
%%%%%%%% Created 01/31/2012  %%%%%%%%%
%Interface with the tracking programs
%Calculate and return migration parameters U,Vx,Vy,P,Px,Py
%Generate a polar plot and histogram of U,Vx,Vy
%Creat an output file,a text file with a name of fname
%Require two input parameters: data and fname
%data is track infomation obtained from the trackikng code(>>l.cellArray)
%fname takes string vector and will be the name of output file 

n= length(data);

%Initialize variables & vectors
U=[];xvel=[]; yvel=[];Px=[]; Py=[];P=[];
stdx=[]; stdy=[];Vx=[]; Vy=[];
i=0; mint=inf;

%Conversion factor, convert pixel to micron
pixpermic=1;

%Time interval between 2 consecutive time frames
timeStep=8;

%Criteria for motile cells
x_cutoff = 2;
y_cutoff =2.5;
% stdx and stdy are used as criteria in this case
% any cells with stdx and stdy smaller than the x_cutoff and y_cutoff
% will be not considered in the analysis, but still on the polar plots


% Creat a text file for recording migration parameters of an individaul cell
% and also the average values over the cell population of the same data set
fid = fopen(['result' fname '.txt'],'w');
fprintf(fid,'\n x= %2.1f  y= %2.1f\n\n',x_cutoff,y_cutoff);
fprintf(fid,'ID  Speed        x-velocity       y-velocity       x-Plength       y-Plength      Plength   Track Duration  \n\n');

%Color vector for polar plot
colrVec = colorInterp(n);

%Loop over entire data set to calculate for migration parameters of each
%cell
for ind=1:n
    
    x = data(ind).X(:,1).*pixpermic;
    y = data(ind).X(:,2).*pixpermic;
    t = data(ind).T.*timeStep;
    
    %Calculate speed, std, persistence length, and velocity
    allSpd(ind) = calSpeed(x,y,t);
    [stdx(ind),stdy(ind)]= calStd(x,y);
    [xPlength(ind),yPlength(ind),Plength(ind)]= calPlength(x,y,t);
    [xvel(ind),yvel(ind)] = calVelocity(x,y,t);
    
    
    tl= t(length(t))-t(1);
    %Screen out non-motile cells
    if stdx(ind) >= x_cutoff && stdy(ind)>= y_cutoff && tl>= 10*timeStep
        i=i+1;
        U = [U allSpd(ind)];
        Vx= [Vx xvel(ind)];
        Vy= [Vy yvel(ind)];
        Px= [Px xPlength(ind)];
        Py= [Py yPlength(ind)];
        P= [P Plength(ind)];
        
        
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
        
        %pause
    end
    
    
    %Figure1: Trajectory plot/ Polar plot
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


fprintf(fid,'# total cells:      %d\n', n);
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

%%%%%%%%%%%%%%%%%%%% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%% calculating for persistence length %%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%% calculating for x,y standard deviation%%%%%%%%%%%%%%%
function [stdx,stdy]= calStd(x,y)

stdx= std(x);
stdy= std(y);

