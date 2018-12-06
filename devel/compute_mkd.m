% computes M^(-1/2) K M(-1/2) d
function f = compute_mkd(d,m,H,Ht,a)

ngll = size(H,1);
displ = reshape(d, ngll,ngll);

displ = displ./sqrt(m);

dU_dxi = Ht * displ;
dU_deta = displ * H;

tmp = a(:,:,1).*dU_dxi + a(:,:,3).*dU_deta;
f = H * tmp;

tmp = a(:,:,3).*dU_dxi + a(:,:,2).*dU_deta;
f = f + tmp * Ht;

f = f./sqrt(m);

f = reshape(f,ngll*ngll,1);
