function aviobj = makemovie(filename, frameset, plot_func, varargin)
% MAKEMOVIE  Create an AVI movie.
%  aviobj = makemovie(filename, frameset, plot_func [, ...])
%
%  Creates an AVI movie file named 'filename.avi' in the current directory.
%  Frame k in the movie is obtained from the contents of the
%  current figure after calling plot_func(frameset{k}).
%  The optional extra parameters are passed directly to avifile.
%
%  Returns the closed avi object handle.
%
%  Example: makemovie('test', states, @(x) plots.tomography(x))

% James D. Whitfield 2009
% Ville Bergholm 2009


if (nargin < 3)
  error('Usage: makemovie(filename, frameset, plot_func, ...')
end

% create an AVI object
aviobj = avifile(filename, varargin{:});

fig = figure('Visible', 'off');

for k = 1:length(frameset)
  plot_func(frameset{k});

  aviobj = addframe(aviobj, fig);

%  F = getframe(fig);   
%  aviobj = addframe(aviobj, F);
end

close(fig)
aviobj = close(aviobj);
