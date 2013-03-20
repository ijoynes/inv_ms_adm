function g_proj = projgr(x, g, lb, ub, nbd)

assert(nargin == 5);
assert(nargout == 1);

assert(isvector(x));
assert(isnumeric(x));
n = length(x);

assert(isvector(lb));
assert(isnumeric(lb));
assert(length(lb)== n);


assert(isvector(ub));
assert(isnumeric(ub));
assert(length(ub)== n);

assert(isvector(nbd));
assert(isnumeric(nbd));
assert(length(nbd)== n);

for i = 1 : n
  assert( 0 <= nbd(i) & nbd(i) <= 3);
end

assert(isvector(g));
assert(isnumeric(g));
assert(length(g)== n);

g_proj = g;
for i = 1 : n
   if nbd(i) ~= 0
       if g(i) < 0
           if nbd(i) >= 2
               g_proj(i) = max((x(i)-ub(i)),g(i));
           end
       else
           if nbd(i) <= 2
               g_proj(i) = min((x(i)-lb(i)),g(i));
           end
       end
   end
end