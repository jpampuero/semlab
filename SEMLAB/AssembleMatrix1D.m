% Assemble a global matrix from its local contributions

function K = AssembleMatrix1D(Kloc,iglob)

[NGLL,NEL] = size(iglob);
nglob = NGLL*NEL;
K = zeros(nglob,nglob);
for e=1:NEL,
  ig = iglob(:,e);
  K(ig,ig) = K(ig,ig) + Kloc(:,:,e);
end 

