function [q] = PunktPolozenie(Q,nr_czlonu,S_A)
q=Q(3*nr_czlonu-2:3*nr_czlonu-1,:);
for i=1:length(q)
    q(:,i)=q(:,i)+Rot(Q(3*nr_czlonu,i))*S_A;
end
end