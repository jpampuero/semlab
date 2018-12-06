function pel = GetLocalValue(prop,e,x)

if isstruct(prop) % a function
  pel = feval(prop.fun,prop.data,e,x);
elseif numel(prop)==1 % a constant
  pel = repmat(prop,size(x));
elseif any(size(prop)==1) % element-wise constant 
  pel = repmat(prop(e),size(x));
else % a table
  pel = prop(:,e);
end

