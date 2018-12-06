% stability of deformed Q4 elements
Ps = [2:8];
Ps = 2;

%x3 = [0.2 0.4 0.6 0.7 0.8 0.9 1:0.2:1.8 2:10]; 
x3 = [0.6 0.7 0.8 0.9 1:0.2:1.8 2:10]; 
y3 = 1.*x3;

coorg = [ [0;0] [1;0] [1;1] [0;1]]; 
omega = zeros(length(x3),length(Ps));
omega_dx = zeros(length(x3),length(Ps));
omega_dx2 = zeros(length(x3),length(Ps));
omega_dx3 = zeros(length(x3),length(Ps));

opts.issym = 1;
opts.isreal = 1;
opts.disp = 0;

for p = 1:length(Ps),

  ngll = Ps(p)+1;
  [xgll,wgll,H] = GetGLL(ngll);
  w2 = wgll * wgll' ;
  Ht = transpose(H);
  a = zeros(ngll,ngll,3);

for k=1:length(x3),
  
  coorg(1,3) = x3(k); 
%  coorg(2,3) = y3(k);
  coorg(1,2) = x3(k); 

  coord = zeros(2,ngll,ngll);
  M = zeros(ngll,ngll);
  DxiDx = zeros(ngll,ngll);
  DxiDz = zeros(ngll,ngll);
  DetaDx = zeros(ngll,ngll);
  DetaDz = zeros(ngll,ngll);
  dxmin =inf;
  dxmin2 =inf;
  dxmin3 =inf;

  for j = 1:ngll,
  for i = 1:ngll,

    shape = q4_elem('shape',xgll(i),xgll(j));
    coord(:,i,j) = coorg * shape';
    if i>1 & j>1

      % based on edges
      dxmin = min( [dxmin, ...
                    norm(coord(:,i,j)-coord(:,i,j-1)), ...
                    norm(coord(:,i,j)-coord(:,i-1,j))] );

      % based on edges and diagonals
      dxmin2 = min( [dxmin2, ...
                    norm(coord(:,i,j)-coord(:,i,j-1)), ...
                    norm(coord(:,i,j)-coord(:,i-1,j)), ...
                    norm(coord(:,i,j)-coord(:,i-1,j-1))/2, ...
                    norm(coord(:,i-1,j)-coord(:,i,j-1))/2] );

      % based on edges and reciprocal lattice
      % exact for rectangular and parallepipedic elements
      d1 = coord(:,i,j)-coord(:,i-1,j-1);
      d2 = coord(:,i,j-1)-coord(:,i-1,j);
      dxmin3 = min( [dxmin3, ...
                    norm(coord(:,i,j)-coord(:,i,j-1)), ...
                    norm(coord(:,i,j)-coord(:,i-1,j)), ...
                    0.5*abs(d1(2)*d2(1)-d1(1)*d2(2))/norm(d1), ...
                    0.5*abs(d1(2)*d2(1)-d1(1)*d2(2))/norm(d2) ]);
% for rectangles: 1/sqrt(  1/norm(coord(:,i,j)-coord(:,i,j-1))^2 +1/norm(coord(:,i,j)-coord(:,i-1,j))^2 ) 
% see Hughes p.517
    end

    dshape = q4_elem('dershape',xgll(i),xgll(j));
    xjac = coorg * dshape;
    M(i,j) = (xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1))*w2(i,j);
    xjaci = inv(xjac);
    DxiDx(i,j)  = xjaci(1,1);
    DxiDz(i,j)  = xjaci(1,2);
    DetaDx(i,j) = xjaci(2,1);
    DetaDz(i,j) = xjaci(2,2);
  end
  end
  a(:,:,1) = M.*( DxiDx.*DxiDx + DxiDz.*DxiDz );
  a(:,:,2) = M.*( DetaDx.*DetaDx + DetaDz.*DetaDz );
  a(:,:,3) = M.*( DxiDx.*DetaDx + DxiDz.*DetaDz );

  % compute highest frequency
  [v,om2] = eigs(@compute_mkd,ngll*ngll,1,'LM',opts, M,H,Ht,a);
  v = reshape(v, ngll,ngll)./sqrt(M);
  omega(k,p) = sqrt(om2);

  omega_dx(k,p) = sqrt(2)*7/3./dxmin;
  omega_dx2(k,p) = 7/3./dxmin2;
  omega_dx3(k,p) = 7/3./dxmin3;

end
end

subplot(311)
plot(x3,omega./omega_dx, '-o')
ylabel('\Omega / \Omega[\Delta x_{min}]')
subplot(312)
plot(x3,omega./omega_dx2, '-o')
ylabel('\Omega / \Omega[\Delta x_{min 2}]')
subplot(313)
plot(x3,omega./omega_dx3, '-o')
ylabel('\Omega / \Omega[\Delta x_{min 3}]')
xlabel('x3')
