classdef celllib < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Cellarray;
        Colonies;
        Namelib;
        Mask;
    end
    
    properties 
        MAX_X_GAP;
        MAX_T_GAP;
        MAX_CC_GAP;
        MIN_CELLCYCLE;    
        MIN_VALID_T_FRAMES;
    end
    
    methods
        function l = celllib(xgap,tgap,ccgap,cycle,tframe)
            l.Cellarray = [];
            l.Colonies = [];
            l.Namelib = zeros(1,5000);
            l.Mask = [];
            if nargin > 0
                l.MAX_X_GAP = xgap;
                l.MAX_T_GAP = tgap;
                l.MIN_CELLCYCLE = cycle;
                l.MAX_CC_GAP = ccgap;
                l.MIN_VALID_T_FRAMES = tframe;
            end
        end
        
        function savecll(celllib, varname)
            eval( [varname '= celllib.Cellarray;'] );
            save( num2str(varname), varname);
        end
        
        function savecly(celllib, varname)
            if ~isempty(celllib.Colonies) && length(celllib.Mask) ~= length(celllib.Colonies)
                nullmask(celllib);
            end
            selectedcolonies = [];
            indices = find(celllib.Mask == 0);
            n = length(indices);
            for i = 1:n
                selectedcolonies = [selectedcolonies celllib.Colonies(indices(i))];
            end
            eval( [varname '= selectedcolonies;'] );
            save( num2str(varname), varname);
        end
        
        function nameindices = assignname(celllib, n)
            a = find(celllib.Namelib == 0);
            nameindices = zeros(1,n);
            for i = 1:n
                celllib.Namelib(a(i)) = 1;
                nameindices(i) = a(i);
            end
        end
        
        function recyclename(celllib, nameindices)
            celllib.Namelib(nameindices) = 0;
        end
        
        function renameall(celllib)
            n = length(celllib.Cellarray);
            celllib.Namelib = zeros(1,5000);
            nameindices = assignname(celllib, n);
            for i = 1:n
                celllib.Cellarray(i).Name = num2str(nameindices(i));
            end
        end
              
        function  disptrack(celllib)
            n = length(celllib.Cellarray)
            tmax = 0;
            for i = 1:n
                c = celllib.Cellarray(i);
                if c.T(length(c.T)) > tmax
                    tmax = c.T(length(c.T));
                end
            end
            for i = 1:n
            %    plainplot(celllib.Cellarray(i));
                 rainbowplot(celllib.Cellarray(i),colormap(hsv(ceil(tmax))));
            %    pause
            end
%            legend(line([0 0], [1 0],'color',[1 0 0], 'linewidth', 2.5 ), 'time');
            
        end
        
        function dispmovie(celllib, to) %to = 0 for common use
            n = length(celllib.Cellarray)
            tmax = 0;
            for i = 1:n
                c = celllib.Cellarray(i);
                if c.T(length(c.T)) > tmax
                    tmax = c.T(length(c.T));
                end
            end
            hold on;
            for t = to+1:tmax
                set(gca, 'color', [0 0 0]);
                axis equal;
                for i = 1:n
                    c = celllib.Cellarray(i);
                    if t >= c.T(1) && t <= c.T(length(c.T))
                        x1 = getposition(c,t);
                        x2 = getposition(c,t-1);
%                         if mod(t,10) == 1
%                             text(x1(1), x1(2),x1(3),num2str(i));
%                         end
                        x = [x1;x2];
                        line(x(:,1),x(:,2),x(:,3), 'color', [1 1 1]);                        
                    end
                end
                pause(1);
                if mod(t,10) == 0
                    text(0,0,0, num2str(t));
                    pause;
                    figure;
                end
            end
        end
        
        function dispcolonies(celllib)
            if ~isempty(celllib.Colonies) && length(celllib.Mask) ~= length(celllib.Colonies)
                nullmask(celllib);
            end
            indices = find(celllib.Mask == 0);
            n = length(indices);
            timespan = [inf, -inf];
            for i = 1:n
                if celllib.Colonies(i).Timespan(1) < timespan(1)
                    timespan(1) = celllib.Colonies(i).Timespan(1);
                end
                if celllib.Colonies(i).Timespan(2) > timespan(2)
                    timespan(2) = celllib.Colonies(i).Timespan(2);
                end
            end
            for i = 1:n
                colonytrack(celllib.Colonies(indices(i)),colormap(hsv(ceil(timespan(2)-timespan(1))+1)));                
            end
        end
        
        function disptrees(celllib)
            n = length(celllib.Colonies);
            maxdep = 0;
            for i = 1:n
                if celllib.Colonies(i).Depth > maxdep
                    maxdep = celllib.Colonies(i).Depth;
                end
            end
            for i = 1:n
                text(i*2^(maxdep+1), celllib.Colonies(i).Rootcell.T(1)-4, num2str(i),'FontSize',12,'BackgroundColor',[.7 .9 .7]);
                treeplot(celllib.Colonies(i),i*2^(maxdep+1));
            end
        end
        
        function disppedigrees(celllib)
            if ~isempty(celllib.Colonies) && length(celllib.Mask) ~= length(celllib.Colonies)
                nullmask(celllib);
            end
            indices = find(celllib.Mask == 0);
            n = length(indices);
            maxdepth = 0;
            for i = 1:n
                if celllib.Colonies(i).Depth > maxdepth
                    maxdepth = celllib.Colonies(i).Depth;
                end
            end
            for i = 1:n
                text(i*2^(maxdepth+1), celllib.Colonies(indices(i)).Rootcell.T(1)-4, num2str(indices(i)),'FontSize',12,'BackgroundColor',[.7 .9 .7]);
                treeplot(celllib.Colonies(indices(i)),i*2^(maxdepth+1));
            end
        end
        
        function nullmask(celllib)
            celllib.Mask = zeros(size(celllib.Colonies));
        end
        
        function addmask(celllib, filename, filetype)
            if size(celllib.Mask) ~= size(celllib.Colonies)
                celllib.Mask = zeros(size(celllib.Colonies));
            end
            if nargin > 1
                cps = getcrops(filename, filetype);
                for i = 1:length(celllib.Colonies)
                    d = zeros(size(cps));
                    d(:,1) = celllib.Colonies(i).Rootcell.X(1,1) >= cps(:,1);
                    d(:,2) = celllib.Colonies(i).Rootcell.X(1,1) <= cps(:,2);
                    d(:,3) = celllib.Colonies(i).Rootcell.X(1,2) >= cps(:,3);
                    d(:,4) = celllib.Colonies(i).Rootcell.X(1,2) <= cps(:,4);
                    d = ~sum(ones(size(d))-d,2);
                    if ~isempty(find(d ~= 0))
                        celllib.Mask(i) = 0;
                    else
                        celllib.Mask(i) = 1;
                    end
                end
            end
        end
        
        function t = readtrack(celllib, filename)
            a = xlsread(filename);
            X = a(:,1:3);
            T = a(:,7);
            ID = a(:,8);
            i = 1;
            n = length(ID);
            while i <= n
                ids = find(ID == ID(i));
                x = X(ids,:);
                [time, IX] = sort(T(ids));
                position = x(IX,:);
                newcell = stemcell(position,time, num2str(assignname(celllib,1)));
                celllib.Cellarray = [celllib.Cellarray newcell];
                i = max(ids)+1;
            end
            tmax = max(T);
            t = linspace(1,tmax,tmax)';
        end
        
        function b = validtrack(celllib, T, t)
            b = 1;
            if T(1) > t(1) %+ celllib.MAX_T_GAP
                if length(T) < celllib.MIN_VALID_T_FRAMES
                    b = 0;
                end
            end
        end
        
        function filter(celllib, time)         
            if isempty(celllib.Cellarray)
                return;
            end
            i = 1;
            n = length(celllib.Cellarray);
            while i <= n
                j = 1;
                while j < length(celllib.Cellarray(i).T)
                    if celllib.Cellarray(i).T(j+1)-celllib.Cellarray(i).T(j) > celllib.MAX_T_GAP
                        celllib.Cellarray(i).T(j+1)
                        newcell = split(celllib.Cellarray(i),celllib.Cellarray(i).T(j+1),num2str(assignname(celllib,1)));
                        celllib.Cellarray = [celllib.Cellarray newcell];
                    end
                    j = j+1;
                end
                i = i+1;
                n = length(celllib.Cellarray)
            end
            i = n;
            while i > 0
                if ~validtrack(celllib, celllib.Cellarray(i).T, time)
                   recyclename(celllib, str2num(celllib.Cellarray(i).Name));
                   celllib.Cellarray(i) = []; 
                end
                i = i-1;
            end
            n = length(celllib.Cellarray);
            nameid = zeros(1,n);
            for i = 1:n
                nameid(i) = str2num(celllib.Cellarray(i).Name);
            end
            recyclename(celllib,nameid);
            nameid = assignname(celllib,n);
            for i = 1:n
               celllib.Cellarray(i).Name = num2str(nameid(i));
            end
        end
        
        function filter2(celllib, time)
            n = length(celllib.Colonies);
            if isempty(celllib.Colonies)
                return;
            end
            i = n;
            while i > 0
                if celllib.Colonies(i).Count == 1 % && celllib.Colonies(i).Rootcell.T(1) > time(1)              
                   for j = 1:length(celllib.Cellarray)
                       if celllib.Cellarray(j) == celllib.Colonies(i).Rootcell
                           recyclename(celllib, str2num(celllib.Cellarray(j).Name));
                           celllib.Cellarray(j) = [];
                           celllib.Cellarray;%
                           break;
                       end
                   end
                   celllib.Colonies(i) = [];
                   celllib.Colonies;%
                end
                i = i-1;
            end
            m = length(celllib.Cellarray);
            nameid = zeros(1,m);
            for i = 1:m
                nameid(i) = str2num(celllib.Cellarray(i).Name);
            end
            recyclename(celllib,nameid);
            nameid = assignname(celllib,m);
            for i = 1:m
               celllib.Cellarray(i).Name = num2str(nameid(i));
            end
        end
        
        function TrackHungarian (celllib, data, time)
           
            nt = length(time);
            if nt == 0 || nt ~= size(data,3)
                return;
            end
            time = [time; time(nt)];
            u = find(~isnan(data(:,1,1)));
            cellarray = stemcell.package(data(u,:,1), time(1), assignname(celllib,length(u)));
            calxnexts(cellarray,time(2));
            celllib.Cellarray = cellarray;
            
            for i = 2:nt
                a = find(~isnan(data(:,1,i)));
                data2 = data(a,:,i);
                           
                [I, J] = MatchHungarian(celllib, celllib.Cellarray, data2, celllib.MAX_X_GAP, time(i));
                for k = 1:length(J)
                     update(celllib.Cellarray(I(k)), data2(J(k),:), time(i), time(i+1));
                end
                                    
                cJ = celllib.complement(J,size(data2,1));           
                for l = 1:length(cJ)
                    newcell = stemcell(data2(cJ(l),:), time(i), num2str(assignname(celllib,1)));                   
                    celllib.Cellarray = [celllib.Cellarray newcell];
                end
                
                NoiseMonitor(celllib, time, i);
                
                calxnexts(celllib.Cellarray,time(i+1));
 
            end
            renameall(celllib);
        end
        
        function [I, J] = MatchHungarian(celllib, array, pset, tolerance, currenttime)
            I = [];
            J = [];
            if isempty(array) || isempty(pset)
                return;
            end
            m = length(array);
            n = size(pset,1);
            LARGE_NUMBER = 10.^10;
            A = zeros(m,n);
            if length(tolerance) == 1
                tolerance  = tolerance .* ones(m,1);
            end
            for i = 1:m
                a = scalingfactor(celllib, array(i));
                if a == 2
                    a 
                end
                tolerance(i) = tolerance(i)*a;
            end
            for i = 1:m
                for j = 1:n
                    A(i,j) = WeightFun(celllib, array(i), pset(j,:), tolerance(i), currenttime, LARGE_NUMBER);
                end
            end
            [C T] = munkres(A);
            for i = 1:m
                j = find(C(i,:));
                if ~isempty(j) && A(i,j) ~= LARGE_NUMBER
                    I = [I i];
                    J = [J j];
                end
            end
        end
        
        function a = scalingfactor(celllib, cell)
            high = length(cell.T);
            low = max(1,high-10);
            vv = cell.V(low:high,:);
            v = zeros(size(vv,1),1);
            for i = 1:high-low+1
                v(i) = norm(vv(i,:));
            end
            vmax = max(v);
            if vmax >= 0.6*celllib.MAX_X_GAP
                a = 2;
            else
                a = 1;
            end
        end
        
        function w = WeightFun(celllib, cell, pos, tol, ct, LARGE_NUMBER)
            n = length(cell.T);
            if ct <= cell.T(n) + celllib.MAX_T_GAP   
                if norm(cell.X(n,:)-pos) <= tol
                    if cell.T(n) < cell.T(1) + celllib.MIN_VALID_T_FRAMES
                        w = norm(cell.X(n,:)-pos)^2 + 4*celllib.MAX_X_GAP^2;
                    else
                        w = norm(cell.X(n,:)-pos)^2;
                    end
                else
                    w = LARGE_NUMBER;
                end
            else
                w = LARGE_NUMBER;
            end
        end

%         function w = WeightFun(celllib, cell, pos, tol, ct, LARGE_NUMBER)
%             n = length(cell.T);
%             if ct <= cell.T(n) + celllib.MAX_T_GAP
%                 CYCLE_AV = 44;
%                 t = (ct-cell.T(1))/CYCLE_AV;
%                 r = norm(cell.X(n,:)-pos)/celllib.MAX_X_GAP;
%                 if r < 0.5
%                     if t < MIN_VALID_T_FRAMES/CYCLE_AV
%                         w = 
%                     else
%                         
%                     end
%                 else if r <= 1
%                     else
%                 end                
%             else
%                 w = LARGE_NUMBER;
%             end
%         end
        
        function NoiseMonitor(celllib,time,currentindex)
            ct = time(currentindex);
            i = length(celllib.Cellarray);
            while i > 0
                if ct <= celllib.Cellarray(i).T(1) + celllib.MIN_VALID_T_FRAMES
                    count = length(celllib.Cellarray(i).T);
                    den = count/(celllib.MIN_VALID_T_FRAMES+1);
                    low = 1/(celllib.MAX_T_GAP+1);
                    p = (2+low)/3;
                    X = length(celllib.Cellarray(i).T);
                    N = length(time(find(time == celllib.Cellarray(i).T(1)):currentindex));
                    rate = binocdf(X,N,p);
                    if rate < 0.3
                        disp(['delete ', num2str(X), '/',num2str(N),' ', num2str(p), ' ', num2str(rate), ' ', num2str(celllib.Cellarray(i).T(1))]);
                        celllib.Cellarray(i) = [];
                    end
                end
                i = i-1;
            end
        end
                
        function Genealogy (celllib, time)
             nt = length(time);
             for i = 1:nt
                 array = celllib.Cellarray;
                 n = length(array);
                 id1 = [];
                 id2 = [];
                 for j = 1:n
                     if time(i) == array(j).T(1) && isempty(array(j).Ancester)
                         id2 = [id2; j];
                     elseif time(i) >= array(j).T(1) && time(i) <= array(j).T(length(array(j).T))                      
%                         if array(j).T(1) == time(1) 
%                             || ...
%                             time(i)-array(j).T(1)>= celllib.MIN_CELLCYCLE
                            id1 = [id1; j];   
%                         end                
                     end
                 end               
                 [I,J] = MatchColony (celllib, id1, id2, time(i));
              
                 for k = 1:length(I)
                     newcell = split(array(I(k)),time(i),num2str(assignname(celllib,1)));
                     celllib.Cellarray = [celllib.Cellarray newcell];
                     insert(array(I(k)).Colony, newcell, array(I(k)),0);
                     insert(array(I(k)).Colony, array(J(k)), array(I(k)));
                 end
                 cJ = celllib.complement2(J,id2)
                 
                 for l = 1:length(cJ)
                     col = colony('');
                     insert(col,array(cJ(l)),[]);
                     celllib.Colonies = [celllib.Colonies col]; 
                 end
             end             
        end
        
        function [I,J] = MatchColony (celllib, oldid, newid, currenttime)
            array = celllib.Cellarray;
            I = [];
            J = [];
            if isempty(oldid) || isempty(newid)
                return;
            end
            m = length(oldid);
            n = length(newid);
            LARGE_NUMBER = 10.^10;
            A = zeros(m,n);
            
            for i = 1:m
                for j = 1:n
                    A(i,j) = WeightFun2(celllib, array(oldid(i)), array(newid(j)), currenttime, LARGE_NUMBER);
                end
            end
            [C T] = munkres(A);
            for i = 1:m
                j = find(C(i,:));
                if ~isempty(j) && A(i,j) ~= LARGE_NUMBER
                    I = [I oldid(i)];
                    J = [J newid(j)];
                    [currenttime oldid(i) newid(j) A(i,j) oldid(find(A(:,j) < LARGE_NUMBER))']
                end
            end  
            
        end
        
%         function w = WeightFun2(celllib, oldcell, newcell, currenttime, LARGE_NUMBER)
%             x1 = getposition(oldcell, currenttime);
%             x2 = getposition(newcell, currenttime);
%             if norm(x1-x2) <= celllib.MAX_CC_GAP
% %                 if newcell.T(1) == 96
% %                     disp('miao');
% %                 end
%                 dt = currenttime - oldcell.T(1);
%                 tmax = oldcell.T(length(oldcell.T)) - oldcell.T(1);
%                 w = 2 - (1-exp(-(tmax-dt)/10)).*(dt <= tmax).*lognpdf(dt/50,0,0.3)+norm(x1-x2)*1E-5;
%      %          w = 0;
%      %           w = norm(x1-x2);
%             else
%                 w = LARGE_NUMBER;
%             end
%         end
%     end
    
      function w = WeightFun2(celllib, oldcell, newcell, currenttime, LARGE_NUMBER)
            x1 = getposition(oldcell, currenttime);
            x2 = getposition(newcell, currenttime);
            CYCLE_AV = 47;
            t = (currenttime-oldcell.T(1))/CYCLE_AV;
            r = norm(x1-x2)/celllib.MAX_CC_GAP;
            a = 0.5;
            b = 0.7;
            if t < b
                if r < a
                    w = a^2 + 20;
                else
                    if r <= 1
                        w = r^2 + 20;
                    else
                        w = LARGE_NUMBER;
                    end
                end
            else
                if r < a
                    w = a^2/(b-1)^4*(t-1)^4;
                else
                    if r <= 1
                        w = r^2;
                    else
                        w = LARGE_NUMBER;
                    end
                end
            end
        end
    end
    
    methods (Static)
               
        function [e t]= readfiles(filename)
            a = xlsread(filename);
            X = a(:,1:3);
            T = a(:,7);
            tmax = max(T);
            e = nan(200,3,tmax);
            for i = 1: length(T)          
                e(min(find(isnan(e(:,1,T(i))))),:,T(i)) = X(i,:);
            end
            t = linspace(1,tmax,tmax)';
        end
        
        function cI = complement(I,n)
            if n ~= 0
                cI = linspace(1,n,n)';
                cI(I) = [];
            else
                cI = [];
            end
        end
        
        function cI = complement2(I,Iwhole)
            cI = [];
            for i = 1:length(Iwhole)
                j = find(I == Iwhole(i));
                if isempty(j)
                    cI = [cI Iwhole(i)];
                end
            end           
        end
        
        function l = test17
            l = celllib (35, 8, 80, 10, 10);
            [e, t] = celllib.readfiles('RA1_Position.xlsx');
            TrackHungarian(l,e,t);
%           dispmovie(l, 0);
%                filter(l,t);
                disptrack(l);
%            disptrack(l);
%                Genealogy(l,t);
%             
%               filter2(l,t);
%              addmask(l,'HL-60 actual_s1_ss.tif','tiff');
%               dispcolonies(l);
%             disptrees(l);
%               figure(2);
%               disppedigrees(l);
        end
        
        function l = Gifts
            l = celllib (20, 3, 80, 10, 10);
            [e, t] = celllib.readfiles('Calib_D4P2.xls');
            TrackHungarian(l,e,t);
            disptrack(l);
        end

        
    end% methods
    
end%class

