function H = compute_observation_matrix(domain_tri, domain_xy, sensor_xy)

nReceptors = size(sensor_xy,1);
nNodes = size(domain_xy,1);
nNodesPerElement = size(domain_tri,2);

% Find the indicies of the containg triangles and the corresponding
% barry centric coordinates.
B = nan(nReceptors,);
triids = nan(nReceptors,1);
for i = 1 : nTris
    A = [xy(tri(i,:),:)'; ones(1,nNodesPerElement)];
    b = [pt';ones(1,nPts)];
    bcc = A\b;
    for j = 1 : nPts
    if all(0<=bcc(:,j)) && all(bcc(:,j)<1)
        triids(j) = i;
        B(j,:) = bcc(:,j)';
    end
    end 
end

% Compute the sparse observation matrix
index = index = 0;
for i = 1 : nReceptors
  for j = 1 : nNodesPerElement
    index = index + 1;
    Hv(index) = B(i,j);
    row(index) = i;
    col(index) =  tri(triids, j);
  end
end
H = sparse(row, col, Hv, nReceptors, nNodes);