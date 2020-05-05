function [V, dVdp] = rayTraceEdges( meshes )
% rayTraceEdges Ray trace the edges of a mesh.

rangex = [-3, 3];
rangey = [-4, 2];
N = 128;
dx = ((rangex(2)-rangex(1))/(N-1));
dy = ((rangex(2)-rangex(1))/(N-1));
xx = rangex(1):dx:rangex(2);
yy = rangey(1):dy:rangey(2);

maxLayers = 10;
countx = zeros(N,1); % number of depth layers at each fragment
intervalcountx = zeros(N,1); % number of intersections at each fragment
county = zeros(N,1); % number of depth layers at each fragment
intervalcounty = zeros(N,1); % number of intersections at each fragment
% don't really need to keep these if we are ray tracing as we can compute 
% the quantities as we go.
bufferx = repmat(struct('depth',1), N, 1 );
buffery = repmat(struct('depth',1), N, 1 );

for i = 1:N
    x = xx(i);
    y = yy(i);
    for k = 1:size(meshes,2)
        mesh = meshes(k);
        pr = reshape( mesh.p, 2, mesh.N )';
        for j = 1:size(mesh.edges,1)
            e = mesh.edges(j,:);
            t = pr(e(2),:) - pr(e(1),:);
            edgeLengthSquared = t(1)^2+t(2)^2;
            n = [ -t(2), t(1) ]; % normal
            c = -n*pr(e(1),:)';
            % line equation
            % n(1)x+n(2)y+c = 0;
            % ray equation
            % r(t) = (x,0) + t*(0,1) = (x,t)
            % substitute one into the other
            % t = - (n(1)x+c) / n(2)
            depthy = - (n(1)*x+c) / n(2);
            % ray equation
            % r(t) = (0,y) + t*(1,0) = (t,y)
            % substitute one into the other
            % t = - (n(1)x+c) / n(2)
            depthx = - (n(2)*y+c) / n(1);
            % where is y along the edge? (same question for x later)
            alphax = (t * ([x;depthy]-pr(e(1),:)')) / edgeLengthSquared;
            if ( alphax >= 0 && alphax <= 1 )
                if ( countx(i) == maxLayers ) 
                    error('maxLayers exceeded! try increasing');
                end
                % intersecting with edge, add into frame buffer
                countx(i) = countx(i) + 1;
                bufferx( i ).edgeix( countx(i) ) = j;
                bufferx( i ).meshix( countx(i) ) = k;
                bufferx( i ).depth( countx(i) ) = depthy;
                bufferx( i ).alpha( countx(i) ) = alphax;
                bufferx( i ).front( countx(i) ) = pr(e(2),1) > pr(e(1),1);
            end
            alphay = (t * ([depthx,;y]-pr(e(1),:)')) / edgeLengthSquared;
            if ( alphay >= 0 && alphay <= 1 )
                if ( county(i) == maxLayers ) 
                    error('maxLayers exceeded! try increasing');
                end
                % intersecting with edge, add into frame buffer
                county(i) = county(i) + 1;
                buffery( i ).edgeix( county(i) ) = j;
                buffery( i ).meshix( county(i) ) = k;
                buffery( i ).depth( county(i) ) = depthx;
                buffery( i ).alpha( county(i) ) = alphay;
                buffery( i ).front( county(i) ) = pr(e(2),2) < pr(e(1),2);
            end
        end
    end
    % have gone through all meshes and all edges... could process volume
    % here!
end

V = 0;
dVdp = zeros( 1, meshes(1,end).DOFIndexOffset + 2*meshes(1,end).N ); 

% find overlapping intervals
% only need to sort if we see 4 or more surfaces at a pixel
for i = 1:N
    if ( county(i) >= 4 ) 
        [buffery(i).depth,ii] = sort(buffery(i).depth);
        buffery(i).edgeix = buffery(i).edgeix(ii);
        buffery(i).meshix = buffery(i).meshix(ii);
        buffery(i).alpha = buffery(i).alpha(ii);
        buffery(i).front = buffery(i).front(ii);
        % count the number of times we go in and out
        % should assume that buffer(i).front(1) is true
        % can perhaps start by assuming no multiple overlaps!
        incount = 0;
        for j = 1:county(i)
            incount = incount + (buffery(i).front(j)*2-1);
            % assume the overlap goes from j to j+1
            if ( incount > 1 )
                intervalcounty(i) = intervalcounty(i) + 1;
                buffery(i).intervalix( intervalcounty(i) ) = j;
                x1 = buffery(i).depth(j);
                x2 = buffery(i).depth(j+1);
                V = V + dy*(x2-x1);
                
                meshixin = buffery(i).meshix(j);
                edgeix = buffery(i).edgeix(j);
                min = meshes(1,meshixin);
                p1ix = min.edges(edgeix, 1);
                p2ix = min.edges(edgeix, 2);
                % this is the x gradient
                dVdp( min.DOFIndexOffset + p1ix*2-1 ) = dVdp( min.DOFIndexOffset + p1ix*2-1 ) - buffery(i).alpha(j) * dy;
                dVdp( min.DOFIndexOffset + p2ix*2-1 ) = dVdp( min.DOFIndexOffset + p2ix*2-1 ) - (1-buffery(i).alpha(j)) * dy;
                
                meshixout = buffery(i).meshix(j+1);
                edgeix = buffery(i).edgeix(j+1);
                mout = meshes(1,meshixout);
                p1ix = mout.edges(edgeix, 1);
                p2ix = mout.edges(edgeix, 2);
                % this is the x gradient
                dVdp( mout.DOFIndexOffset + p1ix*2-1 ) = dVdp( mout.DOFIndexOffset + p1ix*2-1 ) + buffery(i).alpha(j+1) * dy;
                dVdp( mout.DOFIndexOffset + p2ix*2-1 ) = dVdp( mout.DOFIndexOffset + p2ix*2-1 ) + (1-buffery(i).alpha(j+1)) * dy;
            end
        end
    end

    if ( countx(i) >= 4 ) 
        [bufferx(i).depth,ii] = sort(bufferx(i).depth);
        bufferx(i).edgeix = bufferx(i).edgeix(ii);
        bufferx(i).meshix = bufferx(i).meshix(ii);
        bufferx(i).alpha = bufferx(i).alpha(ii);
        bufferx(i).front = bufferx(i).front(ii);
        % count the number of times we go in and out
        % should assume that buffer(i).front(1) is true
        % can perhaps start by assuming no multiple overlaps!
        incount = 0;
        for j = 1:countx(i)
            incount = incount + (bufferx(i).front(j)*2-1);
            % assume the overlap goes from j to j+1
            if ( incount > 1 )
                intervalcountx(i) = intervalcountx(i) + 1;
                bufferx(i).intervalix( intervalcountx(i) ) = j;
                y1 = bufferx(i).depth(j);
                y2 = bufferx(i).depth(j+1);
                V = V + dx*(y2-y1);
                
                meshixin = bufferx(i).meshix(j);
                edgeix = bufferx(i).edgeix(j);
                min = meshes(1,meshixin);
                p1ix = min.edges(edgeix, 1);
                p2ix = min.edges(edgeix, 2);
                % this is the y gradient
                dVdp( min.DOFIndexOffset + p1ix*2 ) = dVdp( min.DOFIndexOffset + p1ix*2 ) - bufferx(i).alpha(j) * dx;
                dVdp( min.DOFIndexOffset + p2ix*2 ) = dVdp( min.DOFIndexOffset + p2ix*2 ) - (1-bufferx(i).alpha(j)) * dx;
                
                meshixout = bufferx(i).meshix(j+1);
                edgeix = bufferx(i).edgeix(j+1);
                mout = meshes(1,meshixout);
                p1ix = mout.edges(edgeix, 1);
                p2ix = mout.edges(edgeix, 2);
                % this is the y gradient
                dVdp( mout.DOFIndexOffset + p1ix*2 ) = dVdp( mout.DOFIndexOffset + p1ix*2 ) + bufferx(i).alpha(j+1) * dx;
                dVdp( mout.DOFIndexOffset + p2ix*2 ) = dVdp( mout.DOFIndexOffset + p2ix*2 ) + (1-bufferx(i).alpha(j+1)) * dx;
            end
        end
    end
end
    

%figure(2)
%cla;
hold on;
for i = 1:N
    x = xx(i);
    for j = 1:countx(i)
        if ( bufferx(i).front(j) ) 
            plot(x,bufferx(i).depth(j),'b.');
        else
             plot(x,bufferx(i).depth(j),'r.');
        end
    end
    for j = 1:intervalcountx(i)
        ix = bufferx(i).intervalix(j);
        y1 = bufferx(i).depth(ix);
        y2 = bufferx(i).depth(ix+1);
        plot([x,x],[y1,y2],'g-');
    end
end

for i = 1:N
    y = yy(i);
    for j = 1:county(i)
        if ( buffery(i).front(j) ) 
            plot(buffery(i).depth(j),y,'c.');
        else
             plot(buffery(i).depth(j),y,'m.');
        end
    end
    for j = 1:intervalcounty(i)
        ix = buffery(i).intervalix(j);
        x1 = buffery(i).depth(ix);
        x2 = buffery(i).depth(ix+1);
        plot([x1,x2],[y,y],'m-');
    end
end

% now draw gradients!
for i=1:size(meshes,2)
    mesh = meshes(1,i);    
    for j=1:mesh.N
        ix = mesh.DOFIndexOffset + j*2-1;
        iy = ix+1;
        dx = dVdp(ix);
        dy = dVdp(iy);
        if ( dx ~= 0 || dy ~= 0 ) 
            px = mesh.p(j*2-1);
            py = mesh.p(j*2);
            plot( [px,px+dx], [py,py+dy], 'k-' );
        end
    end
end
axis([-3,3,-3,2]);

% How to draw this now?