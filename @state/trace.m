function x = trace(s)
% STATE/TRACE  Trace of the state operator.
%  x = trace(s)
%
%  Returns the trace of the state operator of quantum state s.
%  For a pure state this is equal to the squared norm of the state vector.

% Ville Bergholm 2008

if (size(s.data, 2) == 1)
  x = s.data'*s.data;
else
  x = trace(s.data);
end
