function sensorIndex = placeSensors(Dxy,Sxy)
    n = size(Sxy,1);
    nNodes = size(Dxy,1);
    sensorIndex = zeros(n,1);
    for i = 1 : n
        [v,sensorIndex(i)]=min(sum((Dxy - ones(nNodes,1)*Sxy(i,:)).^2,2));
    end
    
