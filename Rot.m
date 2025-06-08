function R = Rot(fi)
%Macierz rotacji dla przypadku 2D - kąt fi wyrażony w radianach
R = [cos(fi) -sin(fi); sin(fi) cos(fi)]; 
end