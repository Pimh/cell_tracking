classdef stemcell < handle
    %A stem cell
    % To be tested: issame,  package
    properties
        X;
        Xnext;
        Traj;
        V;
        T;
        emean;
        estd;
        Ancester;
        Daughters;
        Name;
        Colony;
    end
    
    methods
        function c = stemcell(x, t, name)  % Think of a way to deal with the nasty constructor!!!!!!!!!!!!!!!!!!!
            if nargin > 0
                c.X = x;
                c.T = t;
                vcal(c);
                n = length(c.T);
                c.Name = name;
                c.Colony = [];
                c.Traj = [];
                trajcal(c);
                c.emean = geterrormean(c);
                c.estd = geterrorstd(c);
                c.Xnext = [0 0 0];
            end
        end
        
        function vcal(c)
            if isempty(c.X)
                return;
            end
            c.V(1,:) = [0 0 0];
            n = length(c.T);
            for i=1:n-1
                c.V(i+1,:) = (c.X(i+1,:)-c.X(i,:))./(c.T(i+1)-c.T(i));
            end
        end
        
        function trajcal(c)
            for i = 1:length(c.T)-1
                c.Traj = [c.Traj; stemcell.trajextrapolation(c.X(i,:),c.V(i,:),c.T(i+1)-c.T(i))];
            end
        end
        
        function calxnext(c, tnext)
            n = length(c.T);
            if n == 0
                return;
            else
                c.Xnext = stemcell.trajextrapolation(c.X,c.V,tnext-c.T(length(c.T)));
            end
        end
        
        function calxnexts(c_array, tnext)
            n = length(c_array);
            for i = 1:n
                calxnext(c_array(i),tnext);
            end
        end
        
        function c = update(c,x,t,tnext) %
            n = length(c.T);
            c.T(n+1,1) = t;
            c.X(n+1,:) = x;
            c.V(n+1,:) = (x-c.X(n,:))./(t-c.T(n)); %%%%% to be improve
            c.Traj = [c.Traj; c.Xnext];
            c.emean = geterrormean(c);
            c.estd = geterrorstd(c);
            calxnext(c, tnext);
        end
        
        function newcell = split(oldcell,t0,newcellname)
            indices = find(oldcell.T >= t0);
            newcell = stemcell(oldcell.X(indices,:),oldcell.T(indices),newcellname);
            oldcell.T(indices) = [];
            oldcell.X(indices,:) = [];
            oldcell.V(indices,:) = [];          
        end
        
        function cutoff(c, framenum)
            c.X(1:framenum,:)=[];
            c.Traj(1:framenum-1,:)=[];
            c.V(1:framenum,:)=[];
            c.T(1:framenum,:)=[];
            c.emean = geterrormean(c);  
            c.estd = geterrorstd(c);
        end
        
        function x = getposition(c, t)
            if t >= c.T(1) && t <= c.T(length(c.T))
                id = find( c.T == t );
                if ~isempty(id)
                    x = c.X(id,:);
                else
                    a = max(find(c.T < t));
                    b = min(find(c.T > t));
                    x = [(c.T(b)-t).*c.X(a,:)+(t-c.T(a)).*c.X(b,:)]./(c.T(b)-c.T(a));
                end
            else
                x = [];
            end
        end
        
        function m = geterrormean(c)
            n = length(c.T);
            if n > 1
                dx = c.X(2:n,:)-c.Traj;
                dt = c.T(2:n)-c.T(1:n-1);
                m = mean(dx./[dt dt dt]);
            else
                m = nan;
            end
        end
        
        function Std = geterrorstd(c)
            n = length(c.T);
            if n > 1
                dx = c.X(2:n,:)-c.Traj;
                dt = c.T(2:n)-c.T(1:n-1);
                Std = std(dx./[dt dt dt],0,1);
            else
                Std = nan;
            end
        end
        
        function b = issame(c,x)%x = 1X3
            diff = c.Xnext - x;
            b = (norm(diff)< 1); %Important empiracal parameter!!!!!!!!!
        end
        
        function c2 = copyelement(c1)
            c2 = stemcell(c1.X, c1.T,[c1.Name '_cpy']);
        end
        
        function rainbowplot(c, cmp)  %T minimum >= 1      
            n = length(c.T);
            if n > 0
                %plot3(c.X(1,1), c.X(1,2), c.X(1,3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', cmp(floor(c.T(1)),:), 'MarkerEdgeColor', cmp(floor(c.T(1)),:)); 
                plot3(c.X(1,1)-c.X(1,1), c.X(1,2)-c.X(1,2), c.X(1,3)-c.X(1,3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', cmp(floor(c.T(1)),:), 'MarkerEdgeColor', cmp(floor(c.T(1)),:)); 
                hold on;
                %text(c.X(1,1),c.X(1,2),c.X(1,3),[c.Name ',' num2str(c.T(1)) ',' num2str(c.T(length(c.T)))], 'FontSize', 14);
               % [c.Name ',' num2str(c.T(1)) ',' num2str(c.X(1,1))]
               fprintf('%d , %d \n',c.T(1),c.T(length(c.T)));
            end
            for i = 1:n-1
                a = floor(c.T(i+1));
                b = ceil(c.T(i+1));
                if a ~= b
                    tempcmp = cmp(a,:).*(c.T(i+1)-a)+cmp(b,:).*(b-c.T(i+1));
                else
                    tempcmp = cmp(a,:);
                end
                %line(c.X(i:i+1,1),c.X(i:i+1,2),c.X(i:i+1,3),'color',tempcmp, 'linewidth', 2.5);  
                line(c.X(i:i+1,1)-c.X(1,1),c.X(i:i+1,2)-c.X(1,2),c.X(i:i+1,3)-c.X(1,3),'color',tempcmp, 'linewidth', 2.5);  
            end
            grid on;
            axis equal;
            
        end
        
        function plainplot(c)  %T minimum >= 1  
            figure(1);
            set(gca, 'Color', [0 0 0]);
            n = length(c.T);
            if n > 0
                plot3(c.X(1,1), c.X(1,2), c.X(1,3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1]); 
                hold on;
%                 text(c.X(1,1),c.X(1,2),c.X(1,3),[c.Name ',' num2str(c.T(1)) ',' num2str(c.T(length(c.T)))], 'FontSize', 14);
               % [c.Name ',' num2str(c.T(1)) ',' num2str(c.X(1,1))]
              
            end
            for i = 1:n-1
                 line(c.X(i:i+1,1),c.X(i:i+1,2),c.X(i:i+1,3), 'Color', [1 1 1]);
            end
            grid on;
            axis equal;
            
        end
 
    end
    
    methods(Static) % static methods
                
        function xnext = trajextrapolation(x,v,delta_t)
            s = size(x);
            n = s(1);
          %  xnext = x(n,:) + v(n,:).*delta_t;
            xnext = x(n,:);
        end
        
        function b = isequal(c1,c2)
            b = strcmp(c1.Name, c2.Name);
        end
        
        function c = randomcell(index)
            tmax = 100;
            x = zeros(tmax,3);
            t = (linspace(1,tmax,tmax))';
            x(1,1) = 5*(rand-0.5);
            x(1,2) = 5*(rand-0.5);
            x(1,3) = 5*(rand-0.5);
            for i = 1:tmax-1
                x(i+1,1) = x(i,1) + rand-0.5;
                x(i+1,2) = x(i,2) + rand-0.5;
                x(i+1,3) = x(i,3) + rand-0.5;
            end
            c = stemcell(x,t,['randomcell' num2str(index)]);
            
        %    c.X(tmax,1) = 0;
        %    c.Xnext = [0 0 0];
        %    c.V(tmax,:) = [0 0 0];
        %    c.Xnextdl(1).Data = 0;
        %    c.Xnextdl(2).Data = 0;
        %    c.Xnextdl(3).Data = 0;
        end
        
        function c = randomcell2(index, time)
            n = length(time);
            x = zeros(n,3);
            x(1,1) = 5*(rand-0.5);
            x(1,2) = 5*(rand-0.5);
            x(1,3) = 5*(rand-0.5);
            for i = 1:n-1
                x(i+1,1) = x(i,1) + (rand-0.5).*(time(i+1)-time(i));
                x(i+1,2) = x(i,2) + (rand-0.5).*(time(i+1)-time(i));
                x(i+1,3) = x(i,3) + (rand-0.5).*(time(i+1)-time(i));
            end
            c = stemcell(x,time,['randomcell' num2str(index)]);
        end
        
        function cellarray = package(data, time, nameindices)
            s = size(data);
            n = s(1);
            for i = 1:n
                cellarray(i) = stemcell(data(i,:),time,num2str(nameindices(i))); % Name, velocity....
            end
        end
        
        function test
            c = stemcell([1 -1 0],1,'cell1');
            calxnext(c,2);
            c
            c.X
            c.T
            c.V
            c2 = stemcell([0 0 0; 1 0 0],[1;2],'cell2');
            calxnext(c2,3);
            c2

            c3 = stemcell([-0.2 1.1 3.4; 1 2 3; 4 5 6],[1;2;3],'cell3');
            calxnext(c3,4);
            c3

            c3.X
            c3.T
            c3.V
            update(c3,[10 12 11], 5, 6);
            c3

            c3.X
            c3.T
            c3.V
            c4 = copyelement(c3)
            rainbowplot(c4,colormap(hsv(100)));
            c5 = stemcell.randomcell2(5,[1 2 3 4 5 7 9 20.4 30.7 40.5 60.4 80.1 90 95.4]');
            calxnext(c5,97.4);
            c5

            c5.X
            c5.T
            c5.V
            rainbowplot(c5,colormap(hsv(100)));
            e = celllib.randomSample(3, 20)
            a = stemcell.package(e(:,:,1), 1, [6 7 8])
            calxnexts(a,1);
            a(1)

            a(2)
  
            a(3)

        
        end
        
        function test2
            c3 = stemcell([-0.2 1.1 3.4; 1 2 3; 4 5 6],[1;2;3],'cell3');
            calxnext(c3,4);
            c3

            c3.X
            c3.T
            c3.V
            update(c3,[10 12 11], 5, 6);

            c3.X
            c3.T
            c3.V
        end
        
        function c = test3
            c = stemcell([-0.2 1.1 3.4; 1 2 3; 4 5 6],[1;2;3],'cell3');
            calxnext(c, 5);
            c.Xnext
            c.emean
            c.estd
            geterrormean(c)
            geterrorstd(c)
            update(c,[10 12 11], 5, 6);
            c.emean
            c.estd
            geterrormean(c)
            geterrorstd(c)
        end
    end
    
end

