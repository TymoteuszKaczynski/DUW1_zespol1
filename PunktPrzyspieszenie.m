function [d2q] = PunktPrzyspieszenie(Q,dQ,d2Q,nr_czlonu,S_A)
omega = [0 -1; 1 0];
d2q=d2Q(3*nr_czlonu-2:3*nr_czlonu-1,:);
for i=1:length(d2q)
    d2q(:,i)=d2q(:,i)+(omega*Rot(Q(3*nr_czlonu,i))*S_A*d2Q(3*nr_czlonu,i))-(Rot(Q(3*nr_czlonu,i))*S_A*dQ(3*nr_czlonu,i)^2);
end
end