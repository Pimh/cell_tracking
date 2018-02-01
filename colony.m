classdef colony < handle
    %A colony of stemcells comprises cells from a same ancester.

    
    properties
        Stemcell_array; % A arrary of stemcells
        Rootcell; % The oldest cell in the current colony, and its precedents.
        Depth;
        Count;
        Timespan;
        Name;
        CharaVector;
    end
    
    methods
        function C = colony(name)
            C.Name = name;
            C.Stemcell_array = stemcell;
            C.Rootcell = [];
            C.Depth = 0;
            C.Count = 0;
            C.Timespan = [inf -inf];
            C.CharaVector = [];
        end
        
        function renewtimespan(Colony)
            %Renew the time span of the colony.
            Colony.Timespan = [inf, -inf];
            for i = 1:length(Colony.Stemcell_array)
                if Colony.Stemcell_array(i).T(1) < Colony.Timespan(1)
                     Colony.Timespan(1) = Colony.Stemcell_array(i).T(1);
                end
                if Colony.Stemcell_array(i).T(length(Colony.Stemcell_array(i).T)) > Colony.Timespan(2)
                     Colony.Timespan(2) = Colony.Stemcell_array(i).T(length(Colony.Stemcell_array(i).T));
                end
            end
        end
        
        function b = isInColony(Colony, c) 
            %return whether a cell is a member of the colony.
            n = length(Colony.Stemcell_array);
            for i = 1:n
                if stemcell.isequal(Colony.Stemcell_array(i),c)
                    b = i;
                    return;
                end
            end
            b = 0;
        end
        
        function A = timeview(Colony, num)
            n = length(Colony.Stemcell_array);
            A = zeros(num+1,n);
            for i=1:n
                A(1,i) = str2num(Colony.Stemcell_array(i).Name);
                A(Colony.Stemcell_array(i).T,i) = 1;
            end
        end
        
        function height = getheight(Colony, cell)
            %return the height of the subtree rooted at cell
            if isempty(cell)
                height = 0;
                return;
            end
            n = length(cell.Daughters);
            if n > 0
                daughterheight = zeros(1,n);
                for i = 1:n
                    if ~isempty(cell.Daughters(i))
                        daughterheight(1,i) = getheight(Colony, cell.Daughters(i));
                    end
                end
                height = max(daughterheight)+1;
            else
                height = 1;
            end
        end
        
        function depth = getdepth(Colony, cell)
            % return the depth of the cell in the Colony pedigree tree.
            temp = cell;
            depth = 1;
            while ~isequal(temp,Colony.Rootcell)
                if isempty(temp)
                    depth = 0;
                    return;
                end
                depth = depth+1;
                temp = temp.Ancester;
            end
        end
        
        function insert(Colony, newcell, oldcell, flag)
            %insert a newcell as the descendent of the oldcell in Colony.
            if nargin < 4
                flag = 1;
            end
            if isempty(oldcell)
                Colony.Rootcell = newcell;
                Colony.Count = 1;
                Colony.Depth = 1;
                Colony.Stemcell_array = newcell;
                Colony.Timespan(1) = newcell.T(1);
                Colony.Timespan(2) = newcell.T(length(newcell.T));
                newcell.Colony = Colony;
                Colony.CharaVector = [Colony.CharaVector newcell.T(1)];
                return;
            end
            if ~isInColony(Colony, oldcell)|| isInColony(Colony,newcell)
                return;
            end
            Colony.Stemcell_array = [Colony.Stemcell_array newcell];
            Colony.Count = Colony.Count+1;
            if getdepth(Colony, oldcell) >= Colony.Depth
                Colony.Depth = Colony.Depth+1;
            end
            oldcell.Daughters = [oldcell.Daughters newcell];
            newcell.Ancester = oldcell;
            newcell.Colony = Colony;
            if newcell.T(length(newcell.T)) > Colony.Timespan(2)
                Colony.Timespan(2) = newcell.T(length(newcell.T));
            end
            if flag == 1
                Colony.CharaVector = [Colony.CharaVector newcell.T(1)];
            end
        end
        
        function c2 = copy(c1)
            %return a copy colony
            c2 = colony;
            n = length(c1.Stemcell_array);
            for i = 1:n
                c2.Stemcell_array(i) = copyelement(c1.Stemcell_array(i));
                if stemcell.isequal(c1.Stemcell_array(i),c1.Rootcell)
                    c2.Rootcell = c2.Stemcell_array(i);
                end
            end
            t = 0;
            t = c1.Depth;
            c2.Depth = t;
            v = 0;
            v = c1.Count;
            c2.Count = v;
            w = 0;
            w = c1.Timespan;
            c2.Timespan = w;
        end
        
        function colonytrack(Colony, cmp, cell)
            %plot all the tracks of a colony
            if nargin == 2
                hold off;
                c = Colony.Rootcell;
            end
            if nargin == 3
                c = cell;
            end
            if ~isempty(c)
                hold on;
                rainbowplot(c,cmp);
                if c == Colony.Rootcell
                    text(Colony.Rootcell.X(1,1), Colony.Rootcell.X(1,2), num2str(Colony.Rootcell.Name),'FontSize',12,'BackgroundColor',[.7 .9 .7]);
                end
                for i = 1:length(c.Daughters) 
                    colonytrack(Colony, cmp, c.Daughters(i));
                end
            end
        end
        
        function treeplot(Colony, rootposition, cell, depth, timelow)
            %plot the binary pedigree tree of Colony with the root position
            %at rootposition.
            if nargin == 2
                hold off;
                c = Colony.Rootcell;
                d = Colony.Depth;
                tlow = Colony.Rootcell.T(1);
            else
                c = cell;
                d = depth;
                tlow = timelow;
            end
            rp = rootposition;
            if d == 0
                return;
            end           
            text(rp, tlow, c.Name,'FontSize',12);
            hold on;
            t1 = nan;
            t2 = nan;
            l = length(c.Daughters);
            if l >= 1            
                t1 = c.Daughters(1).T(1);
                if l == 2
                    t2 = c.Daughters(2).T(1);
                end
            end
            thigh = min(t1,t2);
            if l >= 1
                treeplot(Colony,rp-2^(d-1),c.Daughters(1),d-1,thigh); 
                if l == 2
                    treeplot(Colony,rp+2^(d-1),c.Daughters(2),d-1,thigh); 
                end
            end
            if l == 0
                line([rp rp],[tlow c.T(length(c.T))]);
            else
                line([rp rp],[tlow thigh]);
                line([rp-2^(d-1) rp+2^(d-1)],[thigh thigh]);
            end
            plot(rp.*ones(size(c.T)),c.T,'.');
        end

    end%methods
    
    methods(Static)
        function C = randomColony
            C = colony('randomColony');
            c1 = stemcell.randomcell2(1, linspace(1,10,10)');
            c2 = stemcell.randomcell2(2, linspace(10,23,13)');
            c3 = stemcell.randomcell2(3, linspace(10,34,24)');
            c4 = stemcell.randomcell2(4, linspace(23,44,21)');
            c5 = stemcell.randomcell2(5, linspace(23,41,18)');
            c6 = stemcell.randomcell2(6, linspace(41,56,14)');
            c7 = stemcell.randomcell2(7, linspace(34,77,44)');
            c8 = stemcell.randomcell2(8, linspace(34,36,3)');
            insert(C,c1,[]);
            insert(C,c2,c1);
            insert(C,c3,c1);
            insert(C,c4,c2);  
            insert(C,c5,c2);
            insert(C,c6,c5);
            insert(C,c7,c3);
            insert(C,c8,c3);
        end
        
        function test
            C = colony.randomColony
        %    treeplot(C,0);
         %   colonytrack(C, colormap(hsv(100)));
            treeplot(C, 0);
            getheight(C,C.Rootcell)
            getdepth(C,C.Rootcell)
        end
    end
    
end

