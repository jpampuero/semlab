% optimize symplectic scheme with n=4

% lambda --> x
x=[-1:0.0001:1]';
% eq 16 of Omelyan et al (2002)
k= find(2*x.*(-1+8*x-24*x.^2+24*x.^3)>=0);
x=x(k);
% can be + or - in front of sqrt:
chi = (4*x-8*x.^2 - sqrt(2*x.*(-1+8*x-24*x.^2+24*x.^3))) ...
      ./ (12*x-96*x.^3+96*x.^4 );
xi = 1/6-4*x.^2 .*chi;

% our gamma_2 is = their -gamma_5
g2 = -( 1/120-x.*chi.*(1/6-x/2.*(1/6+chi/2-chi.^2-chi.*xi-xi+xi.^2)-chi/6+chi.*xi/2-5*xi/6+xi.^2) ...
     -xi/16 +7*xi.^2/48 -xi.^3/8) ;

% our gamma_4 is = their -gamma_2
g4 = -(1/480-x.*chi/2 .*(1/12-x/6+x.^2.*chi+x.*xi-chi/12-xi/6)-xi.^2/24);

plot(x,g4-g2, 'o')
axis([-inf inf -0.01 0.01])
