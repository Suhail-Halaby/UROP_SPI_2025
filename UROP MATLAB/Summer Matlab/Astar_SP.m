function [path] = Astar_SP(domain_grid,start_pose,end_pose,evaluated,heuristic,plot_limits)
% This function performes Astar between two points, using an evaluated and
% hueristic cost.

% Open points --> A running record of available options
% Closed Points --> A record of our decisions

% finding indicies of start and endpoints

[~, start_xidx] = min(abs(domain_grid(1,:,1) - start_pose(1)));
[~, end_xidx]   = min(abs(domain_grid(1,:,1) - end_pose(1)));

% For Y index
[~, start_yidx] = min(abs(domain_grid(:,1,2) - start_pose(2)));
[~, end_yidx]   = min(abs(domain_grid(:,1,2) - end_pose(2)));

axis equal
axis([ plot_limits(1) plot_limits(2) plot_limits(3) plot_limits(4)] )
%axis([51.9 52.1 -0.1 0.1])
hold on
% point data structure: [X-coord ; Y-coord ; G-cost ; H-cost ; F_cost ; Parent Index]

open = [start_xidx;start_yidx;0;heuristic(end_pose,start_pose);heuristic(end_pose,start_pose);-1];
closed = [];


% Steps in the loop:
% 1. Pick open node with the lowest total cost
% 2. Check if this node is the goal
        % Recover path
% 3. If not, add this point to the list of closed points
% 4. From there, evaluate the new list of open points, ie, the newly added neighbohrs

    
    while ~isempty(open)
    
        % Step 1:
        [~, idx] = min(open(5,:));
        current = open(:,idx);
        
        % 1: finds index that should be used for the current point
    
        % Step 2:
        if isequal(current(1:2), [end_xidx;end_yidx])
            path = recover_path(closed,current);
            return;
        end

        % 2: recovers the final path if the goal has been acheived


        % Step 3:
        closed = [closed,current];
        
        xcl = domain_grid(1,closed(1,:),1);
        ycl = domain_grid(closed(2,:),1,2);
        plot(xcl,ycl,'.b',MarkerSize=10)

        open(:,idx) = [];

        % 3: adds the current list to the set of closed (visited points),
        % wipes the current point from the list of open points

        % Step 4:
        next_opns = check_options(current,domain_grid);

        for n = 1:size(next_opns,2)
            % eliminating the previous closed points from the option pool
            if any(all(next_opns(:,n) == closed(1:2,:),1))
                continue;
            end
            % evaluating costs
            xc = domain_grid(1,current(1),1);
            yc = domain_grid(current(2),1,2);
            plot(xc,yc,'.g',MarkerSize=10)

            xn = domain_grid(1,next_opns(1,n),1);
            yn = domain_grid(next_opns(2,n),1,2);
            plot(xn,yn,'.r',MarkerSize=10)
            

            
            g_new = current(3) + evaluated([xc;yc],[xn;yn]);
            h_new = heuristic([xc;yc],end_pose);
            f_new = g_new + 1.001*h_new;

            % check if the point is in the open set of points
            % keep the one of lower cost
            % Check if neighbor already in open set
            in_open = find(all(next_opns(:,n) == open(1:2,:), 1), 1); % pick first match
            
            if ~isempty(in_open)
                if g_new < open(3, in_open)
                    % update costs & parent
                    open(3, in_open) = g_new;
                    open(4, in_open) = h_new;
                    open(5, in_open) = f_new;
                    open(6, in_open) = size(closed, 2); % current node is parent
                end
            else
                % add new node
                open = [open [next_opns(:,n); g_new; h_new; f_new; size(closed,2)]];
            end

            % finally update the list of opens
            %open = [open [next_opns(:,n); g_new; h_new; f_new; size(closed,2)]];
            % the parent node that this corresponds to is obtained from the
            % instantaneuos length of the closed array

        end

    
    pause(0.0001)
    end
end
    
function neighbors = check_options(current, domain_grid)
    steps = [-1,1,-1,1,-1,1,0,0; 0,0,-1,1,1,-1,-1,1]; % 8-connectivity
    neighbors = [];
    [nx, ny, ~] = size(domain_grid);

    for k = 1:8
        i = current(1) + steps(1,k);
        j = current(2) + steps(2,k);

        % Bounds check
        if i < 1 || i > nx || j < 1 || j > ny
            continue;
        end

        % Obstacle check
        if domain_grid(j,i,3) == 1
            disp('hit at:')
            disp([i,j])
            continue;
        end

        % Prevent diagonal corner cutting
        if abs(steps(1,k)) == 1 && abs(steps(2,k)) == 1
            if domain_grid(current(1), j, 3) == 1 || ...
               domain_grid(i, current(2), 3) == 1
                continue;
            end
        end

        neighbors = [neighbors [i;j]];
    end
end

function path = recover_path(closed, current)
    path = current(1:2);
    parent = current(6);
    while parent > 0
        current = closed(:,parent);
        path = [current(1:2) path];
        parent = current(6);
    end
end