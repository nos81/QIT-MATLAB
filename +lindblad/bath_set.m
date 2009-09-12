function b = bath_set(b, varargin)
% LINDBLAD/BATH_SET  Change bath properties.
%
%   b = bath_set(b, 'cut_type', 'sharp',  'cut_limit', 5)

% Ville Bergholm 2009


while (length(varargin) >= 2)
  prop = varargin{1};
  val = varargin{2};
  varargin = varargin(3:end);

  switch prop
    case 'cut_limit'
      b.cut_limit = val; % == omega_c/omega0

    case 'cut_type'
      b.cut_type = val;
      
    otherwise
      error('Unknown property "%s"', prop)
  end
end

% update cutoff function (at least Octave uses early binding, so when parameters change we need to redefine it)
switch b.cut_type
  case 'sharp'
    b.cut_func = @(x) (abs(x) <= b.cut_limit); % Heaviside theta cutoff
  case 'exp'
    b.cut_func = @(x) exp(-abs(x)/b.cut_limit); % exponential cutoff
  otherwise
    error('Unknown cutoff type "%s"', b.cut_type)
end
b.dH = [];
b.S  = [];
end
