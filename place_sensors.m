function H = place_sensors(domain_tri, domain_xy, sensor_xy)

nReceptors = size(sensor_xy,1);
nNodes = size(domain_xy,1);
trep = TriRep(tri, xy);

[triids, B] = pointLocation(trep, sensor_xy)
nB = size(B,2);
index = index = 0;
for i = 1 : nReceptors
  for j = 1 : nB
    index = index + 1;
    Hv(index) = B(i,j);
    row(index) = i;
    col(index) =  tri(triids, j);
  end
end
H = sparse(row, col, Hv, nReceptors, nNodes);