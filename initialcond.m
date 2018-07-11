function out = initialcond(S,statenames)
  out = cell2struct(num2cell(S{end}(end,:)),statenames,2);
