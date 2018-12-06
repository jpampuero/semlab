% K = BuildKdd_1D(dd,coor,iglob,alpha);
function K = BuildKdd_1D(dd,coor,iglob,alpha);

if ~isstr(dd), 
  error('SEMLAB:BuildKdd_1D:invalidFirstArgument', ...
        'The first argument (%s) should be a 2-chararcter string',dd)      
end                    
if nargin<4, alpha=1; end

[NGLL,NEL] = size(iglob);
nglob = length(coor);

K = zeros(NGLL,NGLL,NEL);
[xgll,wgll,H] = GetGLL(NGLL);

for e=1:NEL, 
  x = coor(iglob(:,e));
  dx_dxi = 0.5*( x(NGLL)-x(1) );
  alphal = GetLocalValue(alpha,e,x);
  W = alphal.*wgll/dx_dxi;
  W = W(:);

  switch dd 
  
    case '00'
      K(:,:,e) = diag(W);

    case '10',
      %K(:,:,e) = diag(W) * H';
      K(:,:,e) = repmat(W,1,NGLL).* H';

    case '01',
      %K(:,:,e) = H * diag(W);
      K(:,:,e) = H .* repmat(W',NGLL,1);

    case '11',
      %K(:,:,e) = H * diag(W) * H';
      K(:,:,e) = H * ( repmat(W,1,NGLL).* H');

    otherwise,
      error('SEMLAB:BuildKdd_1D:unknownDerivativeType', ...
            'The first argument (%s) should be a valid derivative type (00,10,01,11)',dd)
  end

end
