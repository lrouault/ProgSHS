function [ Tnext ] = DF(a,b,c,dt,dxp,dxm)

% dxp : dx plus => pas d'espace à droite
% dxm : dx moins => pas d'espace à gauche

T=[a;b;c];
k=2;

Tnext = T(k);

Tnext = Tnext + dt * (dxp*T(k-1)+dxm*T(k+1)-(dxm+dxp)*T(k)) / (dxm*(dxp^2)/2 + dxp*(dxm^2)/2);                %conduction dx non uniforme





%Tnext=Tnext-4*dt*(T(k)-T0);                          %convection dx uniforme
%Tnext=Tnext-4*dt*sigma*epsilon*(T(k)^4-T0^4);        %rayonnement dx uniforme


end

