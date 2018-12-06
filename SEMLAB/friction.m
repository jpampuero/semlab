%friction.m	Friction coefficient
% mu = friction(u,f)
%   mu 	friction coefficient
%   u	slip
%   f	friction data
function mu = friction(u,f)

% default = linear slip weakening
if ~isfield(f,'kind'), f.kind=1; end

switch f.kind

 case 1 % Linear slip weakening friction law
  W = (f.MUs-f.MUd) ./ f.Dc;
  mu = max(f.MUs - W.*u, f.MUd);

 case 2 %-- Chambon's non linear law --
% for large slip: strength ~ 1/u^p
  u = u ./(f.p*f.Dc);
  mu = f.MUd +(f.MUs-f.MUd) ./ (1+u).^f.p ; 
% With this definition:
% initial slope = (MUs-MUd)/Dc
% if p=0.4, strength drop at slip=2*Dc : 50% 
%	     			  20   : 80%
%	     			  126  : 90%
%            			  4e4  : 99%
% so if Dc=100e-6 m: 90% drop at slip = 1.2 cm
%                    99% ...          = 4 m

 case 3 %-- Abercrombie and Rice non-linear law --
% for large slip: strength_drop ~ u^p
  u = u ./ f.Dc;
  mu = f.MUs -(f.MUs-f.MUd).* ...
                ( u.*(u<=1) + (f.p-1+u.^f.p)/f.p .*(u>1) );
% With this definition:
% initial slope = (MUs-MUd)/Dc
% p =0.3

 otherwise
  error('friction','unknown kind')

end
