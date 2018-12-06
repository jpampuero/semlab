function [Trup,Tpz] = PlotRuptureFront(V,Veps,D,Dc,X,DT)

[NX,NT] = size(V);
Trup = zeros(NX,1);
Tpz  = zeros(NX,1);

for k=1:NX,

  m = find( V(k,:)>Veps );
  if ~isempty(m)
    m = m(1)-1;
    if m>0
      Trup(k) = m+ (Veps-V(k,m))/(V(k,m+1)-V(k,m));
    else
      Trup(k) = Veps/V(k,m+1);
    end
  else
    Trup(k) = NT;
  end

  m = find( D(k,:)<Dc);
  m = m(end);
  if m<NT
    Tpz(k) = m+ (Dc-D(k,m))/(D(k,m+1)-D(k,m));
  else
    Tpz(k) = NT;
  end

end

Trup = Trup*DT;
Tpz = Tpz*DT;

%plot(X,Tpz,X,Trup)
area( X,Tpz,'FaceColor','b')
hold on; area( X,Trup,'FaceColor','w'); hold off
xlabel('X')
ylabel('T')
