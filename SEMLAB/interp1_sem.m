% INTERP1_SEM interpolates a 1D spectral element field 
%
% db = interp1_sem(x,xi [,iglob])	initialize 
% yi = interp1_sem(y,db [,iglob]) 	interpolate
%
% INPUT		x(NGLL,NEL) the 1D SEM mesh node coordinates
%		xi(:)	locations for interpolation
% 		y	1D SEM field 
%			in local storage y(NGLL,NEL)
%			or in global storage y(:)
%		db	database created by an "initialize" call
%		iglob 	local-to-global index table (1D)
%			required if global storage
%
% OUTPUT	db	interpolation database 
%		yi(:)	field values y interpolated at xi
%		
function out = interp1_sem(x_y, xi_db, iglob)

if nargin<3, iglob=[]; end
if isstruct(xi_db)
  out = do_interp1(x_y, xi_db, iglob);  % interpolate
else
  out = init_interp1(x_y, xi_db, iglob); % initialize
end

%------------ initialize -----------------
function db = init_interp1(x,xi,iglob)

nxi = length(xi);

% if x given in global storage
if nnz(size(x)>1)==1 	
  if isempty(iglob), error('interp1_sem:badInput',...
    'In global storage mode you must give the local-to-global index table'); end
  [NGLL,NEL]=size(iglob);
  x1 = x(iglob(1,:));
  x2 = x(iglob(NGLL,:));

% if local storage
else
  [NGLL,NEL]=size(x);
  x1 = x(1,:);
  x2 = x(NGLL,:);
end

xgll = GetGLL(NGLL);
db.e = zeros(nxi);
db.hxi = zeros(NGLL,nxi);
for k=1:nxi,
  e = find( x1<=xi(k) & xi(k)<=x2 ); % element containing xi
  e = e(1);
  eta = 2*(xi(k)-x1(e))/(x2(e)-x1(e)) -1; % local coordinate of xi in e
  db.e(k) = e;
 % interpolation matrix = values of the NGLL Lagrange polynomials at eta
  db.hxi(:,k) = hgll(eta,xgll);
end


%--------------- interpolate -----------------------
function yi = do_interp1(y,db,iglob)

% if y is a vector, it is given in global storage (#global_node)
if nnz(size(y)>1)==1 	
  for k=1:length(db.e),
    yy(:,k) = y(iglob(:,db.e(k))); 
  end

% if y is a matrix, it is in local storage (#gll_node,#element)
else
  yy = y(:,db.e);
end

yi = sum(db.hxi .* yy)';
