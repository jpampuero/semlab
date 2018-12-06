% PlotSeisTrace(x,t,v)
%
% INPUT		x(nx)	X coordinate of receivers
%		t(nt)	time
%		v(nx,nt) seismograms
%
function PlotSeisTrace(x,t,v)

NT = length(t);
ampli = 0.5*max(abs( v(:) ));
disp(sprintf('Amplitude (trace to trace) = %f',ampli))
if ampli>0
  offset = max(abs(diff(x)));
  ampli = offset/ampli;
end
plot(t, v*ampli +repmat(x,1,NT) );
title('SEM Seismograms')
xlabel('Time')
ylabel('Distance')
