function x = trace(s)
% TRACE  Trace of the lmap.
%  x = trace(s)
%
%  Returns the trace of the lmap s.
%  The trace is only properly defined if s.dim{1} == s.dim{2}.

% Ville Bergholm 2011


if ~isequal(s.dim{1}, s.dim{2})
    error('Trace not defined for non-endomorphisms.')
end
    
x = trace(s.data);
