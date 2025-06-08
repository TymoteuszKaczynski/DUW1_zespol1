function [dq] = PunktPredkosc(Q,dQ,nr_czlonu,S_A)
omega = [0 -1; 1 0];
dq=dQ(3*nr_czlonu-2:3*nr_czlonu-1,:);
for i=1:length(dq)
    dq(:,i)=dq(:,i)+(omega*Rot(Q(3*nr_czlonu,i))*S_A*dQ(3*nr_czlonu,i));
end
end