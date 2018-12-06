if ~exist('DEMO','var'), DEMO = 1; end

switch DEMO,
  case 1,

    H = 1;
    MU1 = 0.25;
    MU2 = 1;
    k = [0.1:0.1:4];

    NGLL = 5;
    Z = [0:H/4:H (H+H/2):H/2:20*H]'; 
    mu = MU1*( Z(1:end-1)<H ) + MU2*( Z(1:end-1)>=H );
    K = sem1d_static_stiffness_SH(k,Z,NGLL,mu);

    % from eq 51 of Ampuero et al (2002):
    Theta = atanh(MU1/MU2);
    Ka = MU1*k./tanh( H*k + Theta );

    subplot(1,4,1)
    mutab = [mu mu]'; mutab = mutab(:);
    ztab = [Z(1:end-1) Z(2:end)]'; ztab=ztab(:);
    plot(mutab,ztab)

    subplot(1,4,[2:4])
    plot(k,K,'o', k,Ka)
    xlabel('k'); ylabel('K'); legend('SEM','analytical');
end 
