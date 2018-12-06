% Analyze and plot results from sem1d_test_*.m

% Overview of the wave field 
% Zoom in the figure to see the reflected waves travelling up and down
%
t = (1:NT)*dt;
figure(1)
clf
PlotSeisTrace(OUTx,t,OUTv);

% The signal at x=0 contains many well separated reflected phases
% The k-th reflected phase has travelled a distance = k*2*L
% so in principle it has a delay k*2*L with respect to the source
% 
tt = 2*L/dt; % round-trip (2*L) travel time
tw = floor(2*Ft0/dt);		% window length (in samples)
Nw = floor( (NT-tw)/tt ) +1; 	% number of reflected phases
dist = [0:Nw-1]*2*L ; 		% travel time 
% extract the reflected phases:
u = zeros(Nw,tw);
for k=1:Nw,
  tp = 1+(k-1)*tt;
  u(k,:) = OUTv(1,tp:tp+tw-1);
end
u(2:end,:) = u(2:end,:)/2; % correct for the reflection coeficient

% Taking the source signal as a reference,
% measure the relative time-delay and amplitude misfit (after time-shift)
% for each reflected phase
%
%uref = u(1,:);
uref = src_timef((1:tw)*dt,'ricker', Ff0,Ft0);
clear err err_ampli err_phase
npad=10; % zero-pad, should be > max expected time error
us = zeros(size(u));
for k=1:Nw,
  err(k) = norm( u(k,:)-uref );
  delay = xcorrshift(uref',[zeros(npad,1);u(k,:)';zeros(npad,1)]); 
  delay = delay +1-(npad+1); % =delay+p1-p2, where p1=1 and p2=npad+1
  us(k,:) = specshift(u(k,:),delay)';
  err_ampli(k) = norm( us(k,:)-uref );
  err_phase(k) = delay*dt;
end 
err= err/norm(uref);
err_ampli = err_ampli/norm(uref);

figure(2)
clf
subplot(321)
PlotSeisTrace([0:Nw-1]'*2*L,(0:tw-1)*dt,us);

subplot(322)
ures = u-repmat(uref,Nw,1);
PlotSeisTrace([0:Nw-1]'*2*L,(0:tw-1)*dt,ures);
title('Misfit')

%dist = dist/ (1/Ff0); % travel time normalized by the fundamental period of the source
subplot(323)
plot(dist,err_ampli,'-o', dist,err,'-x')
%xlabel('Distance / \lambda_0')
xlabel('Travel time')
ylabel('Relative RMS error')
legend('Shifted','Non-shifted',2)

subplot(324)
%plot(dist,err_phase*Ff0,'-o')
%xlabel('Travel time \times Ff0')
%ylabel('Relative phase error \times Ff0')
plot(dist,err_phase,'-o')
xlabel('Travel time')
ylabel('Time delay error')

subplot(325)
plot(dist(2:end),err_ampli(2:end)./dist(2:end),'-o', dist(2:end),err(2:end)./dist(2:end),'-x')
xlabel('Travel time')
ylabel('Relative RMS error / travel time')
legend('Shifted','Non-shifted',2)

subplot(326)
plot(dist(2:end),err_phase(2:end)./dist(2:end),'-o')
xlabel('Travel time')
ylabel('Time delay error / travel time')
