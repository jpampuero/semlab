% PREM model by Victor Tsai (May 3 2012)
%
function [vpv,vph,vsv,vsh,rho] = prem_table_ar(depth,period)
% Input depth in array form (km), single period (s)
vpv = 0*depth;
vph = 0*depth;
vsv = 0*depth;
vsh = 0*depth;
rho = 0*depth;

for i=1:length(depth)
    [vpvt,vpht,vsvt,vsht,rhot]=prem_table(depth(i),period);
    vpv(i) = vpvt;
    vph(i) = vpht;
    vsv(i) = vsvt;
    vsh(i) = vsht;
    rho(i) = rhot;
end

function [vpv1,vph1,vsv1,vsh1,rho] = prem_table(depth,period)
% Input depth in km, period in seconds
% Output vp, vs in km/s; rho in g/cm^3

a=6371;
r = a-depth;
x = r/a;

% Looks up parameters at 1s
if r<0 || r>6371
    error('Depth out of bounds. Exiting.');
elseif r>=0 && r < 1221.5
    rho = 13.0885-8.8381*x^2;
    vp = 11.2622-6.3640*x^2;
    vs = 3.6678-4.4475*x^2;
    Qmu = 84.6;
    Qk = 1327.7;
elseif r>=1221.5 && r<3480
    rho = 12.5815-1.2638*x-3.6426*x^2-5.5281*x^3;
    vp = 11.0487-4.0362*x+4.8023*x^2-13.5732*x^3;
    vs = 0;
    Qmu = inf;
    Qk = 57823;
elseif r>=3480 && r<3630
    rho = 7.9565-6.4761*x+5.5283*x^2-3.0807*x^3;
    vp = 15.3891-5.3181*x+5.5242*x^2-2.5514*x^3;
    vs = 6.9254+1.4672*x-2.0834*x^2+0.9783*x^3;
    Qmu = 312;
    Qk = 57823;
elseif r>=3630 && r<5600
    rho = 7.9565-6.4761*x+5.5283*x^2-3.0807*x^3;
    vp = 24.9520-40.4673*x+51.4832*x^2-26.6419*x^3;
    vs = 11.1671-13.7818*x+17.4575*x^2-9.2777*x^3;
    Qmu = 312;
    Qk = 57823;
elseif r>=5600 && r<5701
    rho = 7.9565-6.4761*x+5.5283*x^2-3.0807*x^3;
    vp = 29.2766-23.6027*x+5.5242*x^2-2.5514*x^3;
    vs = 22.3459-17.2473*x-2.0834*x^2+0.9783*x^3;
    Qmu = 312;
    Qk = 57823;
elseif r>=5701 && r<5771
    rho = 5.3197-1.4836*x;
    vp = 19.0957-9.8672*x;
    vs = 9.9839-4.9324*x;
    Qmu = 143;
    Qk = 57823;
elseif r>=5771 && r<5971
    rho = 11.2494-8.0298*x;
    vp = 39.7027-32.6166*x;
    vs = 22.3512-18.5856*x;
    Qmu = 143;
    Qk = 57823;
elseif r>=5971 && r<6151
    rho = 7.1089-3.8045*x;
    vp = 20.3926-12.2569*x;
    vs = 8.9496-4.4597*x;
    Qmu = 143;
    Qk = 57823;
elseif r>=6151 && r<6291
    rho = 2.6910+0.6924*x;
    vpv=0.8317+7.2180*x;
    vph=3.5908+4.6172*x;
    vsv = 5.8582-1.4678*x;
    vsh = -1.0839+5.7176*x;
    Qmu = 80;
    Qk = 57823;
    eta = 3.3687-2.4778*x;
elseif r>=6291 && r<6346.6
    rho = 2.6910+0.6924*x;
    vpv = 0.8317+7.2180*x;
    vph = 3.5908+4.6172*x;
    vsv = 5.8582-1.478*x;
    vsh = -1.0839+5.7176*x;
    Qmu = 600;
    Qk = 57823;
    eta = 3.3687-2.4778*x;
elseif r>=6346.6 && r<6356
    rho = 2.900;
    vp = 6.800;
    vs = 3.900;
    Qmu = 600;
    Qk = 57823;
elseif r>=6356 && r<6368
    rho = 2.600;
    vp = 5.800;
    vs = 3.200;
    Qmu = 600;
    Qk = 57823;
elseif r>=6368 && r<=6371
    rho = 1.020;
    vp = 1.450;
    vs = 0;
    Qmu = inf;
    Qk = 57823;
else
    error('Invalid input. Exiting');
end

if r<6151 || r>=6346.6
    vpv = vp;
    vph = vp;
    vsv = vs;
    vsh = vs;
    eta = 1;
end

E = 4/3*((vsv+vsh)/(vpv+vph))^2;
vsv1 = vsv*(1-log(period)/(pi*Qmu));
vsh1 = vsh*(1-log(period)/(pi*Qmu));
vpv1 = vpv*(1-log(period)/pi*((1-E)/Qmu+E/Qk));
vph1 = vph*(1-log(period)/pi*((1-E)/Qmu+E/Qk));
    
