function [Rsqr,count]= motility1(fname,refFile,type,colr)
%Read excel file
%Include reference file
%Display info of every single cell that is qualified as a motile cell
%Result is printed on the different files for different device positions

if type==3
    xc=1;yc=2;tc=7;IDc=8;
elseif type ==2
    xc=4;yc=5;tc=3;IDc=2;
end
[cells,nc]=loadCellData(fname,refFile,xc,yc,tc,IDc);


pixpermic=1;

colrVec = colorInterp(nc);
Rsqr=zeros(1,121);
count=zeros(1,121);
for ind=1:nc
    r=[];x=[];y=[];t=[];
    [x,y,t]=getCell(cells,ind);
    x = x.*pixpermic;
    y = y.*pixpermic;
    
    
    r(1)=0;
    for i=2:length(x)
        
            r(i) = (x(i)-x(1))^2+(y(i)-y(1))^2;
            Rsqr(t(i))= Rsqr(t(i))+r(i);
            count(t(i))= count(t(i))+1;
        
    end
    
    
    
%   hold on      
%   figure(1)
%   plot(t,r,'color',colrVec(ind,:));     
    
end

meanRsqr = Rsqr./count;
hold on
figure(1)
plot((1:121),meanRsqr,colr)


% figure(5);plot((1:193)*5,(rs3+rs4)./(t3+t4),'r');hold on;plot((1:192)*5,(rs1+rs2)./(t1+t2),'b');
% title('Motility')
% ylabel('<r^2> (um^2)')
% xlabel('time(min)')


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
fclose(Fcell) ;  