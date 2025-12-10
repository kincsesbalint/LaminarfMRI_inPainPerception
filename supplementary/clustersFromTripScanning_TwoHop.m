function [clustersMarkedOnly, clustersWithBridges] = clustersFromTripScanning_TwoHop(TRIP, ORIGINALVERTICES)
% Cluster marked vertices using 1-hop and 2-hop connectivity derived by scanning TRIP.
% TRIP: T-by-3 triangle vertex indices (1-based).
% ORIGINALVERTICES: vector of marked vertex IDs.
% clustersMarkedOnly: cell array of clusters containing only marked vertices.
% clustersWithBridges: cell array of clusters including unmarked "bridge" vertices.

    % Basic checks
    if isempty(TRIP) || size(TRIP,2) ~= 3
        error('TRIP must be a T-by-3 matrix of vertex indices.');
    end
    marked = unique(ORIGINALVERTICES(:));
    if isempty(marked)
        clustersMarkedOnly = {};
        clustersWithBridges = {};
        return;
    end

    % Number of vertices (assume 1..N)
    n = max(max(TRIP(:)), max(marked));

    % Precompute neighbor lists from TRIP (undirected edges)
    neighbors = buildNeighborsFromTrip(TRIP, n);

    % Map vertex ID -> index in marked list (0 if not marked)
    m = numel(marked);
    posInMarked = zeros(n,1);
    posInMarked(marked) = 1:m;

    visitedMarked = false(m,1);
    clustersMarkedOnly = {};

    % BFS over expanded adjacency (1-hop OR 2-hop) restricted to marked vertices
    for i = 1:m
        if visitedMarked(i), continue; end
        v = marked(i);
        comp = v;
        visitedMarked(i) = true;

        queue = v;
        head = 1;
        while head <= numel(queue)
            u = queue(head);
            head = head + 1;

            % 1-hop neighbors of u
            neigh1 = neighbors{u};
            if isempty(neigh1), neigh1 = []; end

            % Direct marked neighbors
            idx1 = posInMarked(neigh1);
            mask1 = idx1 ~= 0 & ~visitedMarked(idx1);
            toAdd1 = neigh1(mask1);

            % 2-hop candidates: neighbors of neighbors (via any intermediate w)
            if ~isempty(neigh1)
                % Union of all neighbors of neigh1
                neigh2_all = unique(cell2mat(neighbors(neigh1)));
                % Remove self
                neigh2_all(neigh2_all == u) = [];
            else
                neigh2_all = [];
            end
            idx2 = posInMarked(neigh2_all);
            mask2 = idx2 ~= 0 & ~visitedMarked(idx2);
            toAdd2 = neigh2_all(mask2);

            % Add new marked vertices discovered via 1-hop or 2-hop
            toAdd = unique([toAdd1; toAdd2]);
            if ~isempty(toAdd)
                queue = [queue; toAdd]; %#ok<AGROW>
                comp = [comp; toAdd]; %#ok<AGROW>
                visitedMarked(posInMarked(toAdd)) = true;
            end
        end
        clustersMarkedOnly{end+1} = unique(comp); %#ok<AGROW>
    end

    % Build clusters including "joining" vertices (unmarked bridges)
    clustersWithBridges = cell(size(clustersMarkedOnly));
    markedMask = false(n,1); markedMask(marked) = true;

    for k = 1:numel(clustersMarkedOnly)
        compMarked = clustersMarkedOnly{k};

        % Candidate bridge vertices: neighbors of any marked in the component
        if isempty(compMarked)
            clustersWithBridges{k} = compMarked;
            continue;
        end
        compNeighbors = unique(cell2mat(neighbors(compMarked)));

        % Count, for each candidate vertex w, how many marked vertices in the component it touches
        posInCand = zeros(n,1);
        posInCand(compNeighbors) = 1:numel(compNeighbors);
        touchCount = zeros(numel(compNeighbors),1);

        for t = 1:numel(compMarked)
            nei = neighbors{compMarked(t)};
            if isempty(nei), continue; end
            idx = posInCand(nei);
            idx = idx(idx ~= 0);
            touchCount(idx) = touchCount(idx) + 1;
        end

        % Bridge vertices: unmarked and adjacent to >= 2 marked vertices in this component
        isBridge = touchCount >= 2;
        bridgeVerts = compNeighbors(isBridge & ~markedMask(compNeighbors));

        clustersWithBridges{k} = unique([compMarked; bridgeVerts]);
    end
end

% Helper: build neighbor lists by scanning TRIP
function neighbors = buildNeighborsFromTrip(TRIP, n)
    % Edges from triangles
    E = [TRIP(:,[1 2]); TRIP(:,[2 3]); TRIP(:,[3 1])];
    E(E(:,1) == E(:,2), :) = [];     % remove self-edges if any
    undirected = [E; fliplr(E)];     % make undirected pairs
    neighbors = accumarray(undirected(:,1), undirected(:,2), [n,1], @(x){unique(x)});
    % Ensure cells exist for vertices with no neighbors
    for v = 1:n
        if isempty(neighbors{v})
            neighbors{v} = [];
        end
    end
end

