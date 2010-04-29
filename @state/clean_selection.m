function sys = clean_selection(s, sys)
% CLEAN_SELECTION  Internal helper, makes a subsystem set unique and sorted.

  sys = intersect(1:length(s.dim), sys);
