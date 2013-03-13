function s = place_sources(tri, xy, Sxy, Sm)
    nNodes = size(xy,1);
    nTris  = size(tri,1);
    nSources = size(Sxy,1);
    
    sourceIndex = placeSensors(xy,Sxy);
    
    s = zeros(nNodes,1);
    for i = 1 : nSources
        s_temp = zeros(nNodes,1);
        s_temp(sourceIndex(i)) = 1;
        m = 0;
        for j = 1 : nTris
            m = m + det([ones(3,1) xy(tri(j,:),:)])/6*sum(s_temp(tri(j,:)));
        end
        s(sourceIndex(i)) = s(sourceIndex(i)) + Sm(i)/m;
    end