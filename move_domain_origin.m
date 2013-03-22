function x = move_domain_origin(x,origin)
  assert( nargin == 1 || nargin == 2);
  assert( nargout == 1);

  assert(isnumeric(x));
  assert(ismatrix(x));
  assert(all(~isnan(x)));

  assert(isnumeric(origin));
  assert(isvector(origin));
  assert(all(~isnan(origin)));

  assert(size(x,2) == length(origin));

  nDims = size(x,2);

  if nargin == 1              % If no origin was supplied use a default 
    origin = zeros(1,nDims);  %   value of a zero vector for the origin.
  elseif size(origin,1) > 1   % If the origin was supplied as a column 
      origin = origin'        %   vector then change it to a row vector.
  end

  nNodes = size(x,1);
  x = x + (origin - min(x)) * ones(nNodes,1);