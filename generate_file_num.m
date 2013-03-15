function file_num = generate_file_num(n, nt)
%generate_file_num returns a string of digits representing the file 
% number of a time step or iteration data file.  The string of digits
% is sized to accommodate n = nt-1 and padded with leading zeros if 
% less than the allocated string length is required to represent the
% file number.
%
% INPUTS:
%   n   integer representing the current time step or iteration
%   nt  integer representing the total number of time steps or
%         iterations
%
% OUTPUTS:
%   file_num  string of digits representing the file number of a time
%               step or iteration data file.
%
% EXAMPLE:
%   file_num = generate_file_num(42, 100) will return file_num = '42'
%   file_num = generate_file_num(42, 101) will return file_num = '042'

assert( nt >= 1 );  % there must be at least 1 time step or iteration
assert(  n >= 0 );  % the time step or iteration must be positive
assert(  n < nt );  % the time step or iteration must be less than the
                    %   total number of time steps or iterations
if nt == 1
  file_num = '0';
else  % nt > 1
  file_num = repmat('0', 1, ceil(log10(nt)) - ceil(log10(n+1)));
  if n > 0
    file_num = [file_num num2str(n)];
  end
end