% indx = Init2dSnapshot(iglob)
%
% INPUT		iglob(ngll,ngll,nel) local to global numbering map
% OUTPUT	indx(:,4) vertices of each GLL cell
%
function indx = Init2dSnapshot(iglob)

[NGLL,NGLL,NEL] = size(iglob);
indx = zeros(NEL*(NGLL-1)^2, 4);
ip=0;
for e=1:NEL,
  for i=1:NGLL-1,
  for j=1:NGLL-1,
    ip = ip+1;
    indx(ip,:) = [iglob(i,j,e) iglob(i+1,j,e) iglob(i+1,j+1,e) iglob(i,j+1,e)];
  end
  end
end
