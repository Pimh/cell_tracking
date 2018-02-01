function calibratedPos(fname,refname)

R = xlsread(refname);
X= R(:,1);Y=R(:,2);Z=R(:,3);
offset= [X-X(1) Y-Y(1) Z-Z(1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = xlsread(fname);
Pos= D(:,1:3);
T= D(:,7);
ID=1;

[m,n]=size(Pos);
corrected_Pos=zeros(m,n+4);
for t=1:120
    next=0;
    while ~next && ID<=m
        if T(ID)==t
            corrected_Pos(ID,1:3)=Pos(ID,:)-offset(t,:);
            ID=ID+1;
        else
            next=1;
        end
    end
end

corrected_Pos(:,7)= T;
xlswrite(['Calib_' fname],corrected_Pos);