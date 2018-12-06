function [eigval,eigvec,x]=sem1d_two_layers(RefinementRatio,adapted,ngll)

% medium 1, 4*faster: 4*RefinementRatio elements in x=[0:4/5]
% medium 2, 4*slower: 4*RefinementRatio elements in x=[4/5:1]

c1=4;c2=1;
rho1=1;rho2=1;
mu1=rho1*c1^2;mu2=rho2*c2^2;
L=1;H=0.8;

nel1 = RefinementRatio*4;
nel2 = RefinementRatio*4;
nel = nel1+nel2;
if adapted
  Hm=H;
  dxe1 = Hm/nel1;
  dxe2 = (L-H)/nel2;
  mu.fun=@ewise;
  mu.data=[ mu1*ones(nel1,1); mu2*ones(nel2,1)] ;
else
 %material boundary in the middle of an element
  dxe2 = (L-H)/(nel2-0.5);
  Hm = H-0.5*dxe2;
  dxe1 = Hm/nel1;
  mu.fun=@pwise;
  mu.data(1)=mu1; mu.data(2)=mu2; mu.data(3)=H;
end
xmesh = [ (0:nel1)*dxe1, Hm+(1:nel2)*dxe2 ]'; 
rho.fun=@ewise;
rho.data=ones(nel,1);


[eigvec,eigval,x]=sem1d_modal_analysis(xmesh,ngll,'NN',rho,mu);
neig=length(eigval);
k=(0:neig-1)';
if abs(H/c1-(L-H)/c2)<eps
  anaval=pi/2*c1/H*k;
else
  anaval=zeros(neig,1);
  for i=1:neig,
    anaval(i) = fzero(@delta,eigval(i),[],mu1,mu2,c1,c2,H,L);
  end
end

subplot(211)
plot(k,eigval,'o-',k,anaval,'--')
subplot(212)
semilogy(k(2:end),abs(eigval(2:end)-anaval(2:end))./anaval(2:end))

eigval = [eigval,anaval];

%-----------------
% constant inside each element = data(e)
function prop=ewise(data,e,x)
prop=repmat(data(e),length(x),1);

%-----------------
% constant with boundary inside an element
% = data(1) for x<=data(3), else data(2)
function prop=pwise(data,e,x)
xb=data(3);
prop= (x<=xb)*data(1) + (x>xb)*data(2);

%-----------------
function d = delta(w,mu1,mu2,c1,c2,H,L)
the1=w*H/c1;
the2=w*(L-H)/c2;
d= mu1/c1*sin(the1)*cos(the2) + mu2/c2*sin(the2)*cos(the1);
