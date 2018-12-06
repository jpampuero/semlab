%SETTIMESTEP_1D sets a stable timestep for explicit integration 
%
% dt = SetTimestep_1d(CFL, coor,iglob, vs)
%
% INPUT		CFL	Courant stability number 
%		coor(:) Non-redundant list of coordinates of the global GLL nodes
%		iglob	Local-to-global index map (see mesh1d)
%		vs	Wave velocity, in local storage (#igll,#elem)
%
% OUTPUT	dt	time-step = CFL*min(dx/vs)
%		

function dt = SetTimestep_1d(CFL, coor,iglob, vs)

[NGLL,NEL] = size(iglob);
if nargin<4, vs=1; end
if numel(vs)==1, vs=repmat(vs,NGLL,NEL); end

dt = Inf;
for e=1:NEL,
  vloc = 0.5*( vs(1:NGLL-1,e)+ vs(2:NGLL,e) );
  dx = abs(diff( coor(iglob(:,e)) ));
  dt = min(dt, min(dx./vloc));
end
dt = CFL*dt;
