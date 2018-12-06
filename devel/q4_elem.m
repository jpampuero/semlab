function out = q4_elem(mode,s,t)

switch lower(mode)
  case 'shape'

  sp = s + 1;
  sm = s - 1;
  tp = t + 1;
  tm = t - 1;

  out(1) = sm * tm;
  out(2) = - sp * tm;
  out(3) = sp * tp;
  out(4) = - sm * tp;
  
  out = out/4;

  case 'dershape'

  sp = s + 1;
  sm = s - 1;
  tp = t + 1;
  tm = t - 1;

  out(1,1) = tm;
  out(2,1) = - tm;
  out(3,1) =  tp;
  out(4,1) = - tp;

  out(1,2) = sm;
  out(2,2) = - sp;
  out(3,2) =  sp;
  out(4,2) = - sm;

  out = out/4;

end
