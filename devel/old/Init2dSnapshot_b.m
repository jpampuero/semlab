% indx = Init2dSnapshot_b(NGLL)
%
% PURPOSE	Use it with 'local storage scheme'
% INPUT		NGLL	number of GLL nodes per edge
% OUTPUT	indx(:,4) vertices of each GLL cell 
%			in a reference element
%
function indx = Init2dSnapshot_b(NGLL)

iglob = reshape( (1:NGLL*NGLL)', NGLL,NGLL);
indx = zeros((NGLL-1)^2, 4);
ip=0;
for i=1:NGLL-1,
for j=1:NGLL-1,
  ip = ip+1;
  indx(ip,:) = [iglob(i,j) iglob(i+1,j) iglob(i+1,j+1) iglob(i,j+1)];
end
end
