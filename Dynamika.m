
%START.m
fprintf('Krok czasowy:\n')
dT = 0.1; %input("Wpisz dt:\n");
fprintf('Czas symulacji [s]:\n')
T_Max = 30; %input("Wpisz t maksymalne:\n")

[Q, DQ, D2Q] = main(dT,T_Max);


function [Q,dQ,d2Q]=main(dT,T_Max)
czas = 0:dT:T_Max;

[Czlony, ParyObrotowe, ParyPostepowe, WymPost, Masy, Sily, Sprezyto_tlumiace] = WczytajDane();

Y0=[Czlony,Czlony*0];
OPTIONS = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[T,Y]=ode45(@(t,Y) Uklad_ode(t,Y,ParyObrotowe,ParyPostepowe,WymPost,Masy,Sprezyto_tlumiace,Sily),czas,Y0,OPTIONS);
Y=Y';

Q=Y(1:size(Czlony,1),:);
dQ=Y(size(Czlony,1)+1:size(Czlony)*2,:);
d2Q=Q.*0;
for i=2:size(Y,2)-1
    d2Q_temp=Uklad_ode(0,[Q(:,i);dQ(:,i)],ParyObrotowe,ParyPostepowe,WymPost,Masy,Sprezyto_tlumiace,Sily);
    d2Q(:,i+1)=d2Q_temp(size(d2Q_temp)/2+1:size(d2Q_temp));
end
end
function [Fq] = MacierzJacobiego(Q, ParyObrotowe, ParyPostepowe, WymPost)
Omega = [0 -1; 1 0];
Fq = zeros(size(ParyObrotowe,1)+size(ParyPostepowe,1)+size(WymPost,1),length(Q));
indeks_macierz =1;

% Pary obrotowe
for i = 1:size(ParyObrotowe,1)
    for j = 1:2
        if ParyObrotowe(i,j)~=0
            %Wzór 2.31 ze skryptu
            % Fq(1:2, nr.pary obrotowej*3-[2 1])=[1 0;0 1]
            Fq(indeks_macierz:indeks_macierz+1,(-2:-1)+(ParyObrotowe(i,j)*3))=eye(2)*((-1)^(j+1));
            %Wzór 2.32 ze skryptu
            % Fq(1:2, nr.pary obrotowej*3)= Omega*Rot(fi)*sA
            Fq(indeks_macierz:indeks_macierz+1,ParyObrotowe(i,j)*3)=Omega*Rot(Q(3*ParyObrotowe(i,j)))*ParyObrotowe(i,(1:2)+(j*2))'*((-1)^(j+1));
        end
    end
    indeks_macierz=indeks_macierz+2;
end

% Pary postepowe
for i = 1:size(ParyPostepowe,1)
    if ParyPostepowe(i,1)~=0
        ri=Q(3*ParyPostepowe(i,1)-2:3*ParyPostepowe(i,1)-1);
        Roti=Rot(Q(3*ParyPostepowe(i,1)));
    else
        ri=[0;0];
        Roti=Rot(0);
    end
    if ParyPostepowe(i,2)~=0
        rj=Q(3*ParyPostepowe(i,2)-2:3*ParyPostepowe(i,2)-1);
        Rotj=Rot(Q(3*ParyPostepowe(i,2)));
    else
        rj=[0;0];
        Rotj=Rot(0);
    end

    for j = 1:2
        if ParyPostepowe(i,j)~=0
            %Wzór 2.44 ze skryptu
            Fq(indeks_macierz,ParyPostepowe(i,j)*3)=(-1)^(j+1);
            %Wzór 2.49 ze skryptu
            Fq(indeks_macierz+1, ParyPostepowe(i,j)*3-2:ParyPostepowe(i,j)*3-1) = (Rotj * [-ParyPostepowe(i,5);ParyPostepowe(i,4)])'*((-1)^j);
            if j==1
                %Wzór 2.50 ze skryptu
                Fq(indeks_macierz+1, ParyPostepowe(i,j)*3) = -(Rotj * [-ParyPostepowe(i,5);ParyPostepowe(i,4)])' * Omega * Roti * ParyPostepowe(i,6:7)';
            else
                Fq(indeks_macierz+1, ParyPostepowe(i,j)*3) = -(Rotj * [-ParyPostepowe(i,5);ParyPostepowe(i,4)])' * Omega * (rj-ri-Roti *ParyPostepowe(i,6:7)');
            end
        end
    end
    indeks_macierz=indeks_macierz+2;
end

for ii = 1:size(WymPost,1)
    i=WymPost(ii,1);
    if ParyPostepowe(i,1)~=0
        ri=Q(3*ParyPostepowe(i,1)-2:3*ParyPostepowe(i,1)-1);
        Roti=Rot(Q(3*ParyPostepowe(i,1)));
    else
        ri=[0;0];
        Roti=Rot(0);
    end
    if ParyPostepowe(i,2)~=0
        rj=Q(3*ParyPostepowe(i,2)-2:3*ParyPostepowe(i,2)-1);
        Rotj=Rot(Q(3*ParyPostepowe(i,2)));
    else
        rj=[0;0];
        Rotj=Rot(0);
    end

    for j = 1:2
        if ParyPostepowe(i,j)~=0
            %Wzór 2.49 ze skryptu
            Fq(indeks_macierz, ParyPostepowe(i,j)*3-2:ParyPostepowe(i,j)*3-1) = (Rotj * ParyPostepowe(i,4:5)')'*((-1)^j);
            if j==1
                %Wzór 2.50 ze skryptu
                Fq(indeks_macierz, ParyPostepowe(i,j)*3) = -(Rotj * ParyPostepowe(i,4:5)')' * Omega * Roti * ParyPostepowe(i,6:7)';
            else
                Fq(indeks_macierz, ParyPostepowe(i,j)*3) = -(Rotj * ParyPostepowe(i,4:5)')' * Omega * (rj-ri-Roti *ParyPostepowe(i,6:7)');
            end
        end
    end
    indeks_macierz=indeks_macierz+1;
end

% if(det(Fq)==0)
%     error("Błąd w macierzy Jakobiego!")
% end
end
function q=NewRaph(Q, ParyObrotowe, ParyPostepowe, WymPost, t)
% q=NewRaph(q0,t)
%     Rozwiązanie układu równań metodą Newtona-Raphsona
%     
q=Q;
F=ones(length(q),1);
iter=1;
while( (norm(F)>10e-9) & (iter < 25) )
    F=Wiezy(q, ParyObrotowe, ParyPostepowe, WymPost, t);
    Fq=MacierzJacobiego(q, ParyObrotowe, ParyPostepowe, WymPost);
    q=q-Fq\F;
    iter=iter+1;
end
if iter >=25
    error('BŁĄD: Po 25 iteracjach Newtona_Raphsona nie uzyskano zbieżności');
    q=Q;
end
end
function [dQ] = Predkosc(Q, ParyObrotowe, ParyPostepowe,WymPostepowe,czas)

FT = zeros(length(Q), 1);
indeks_macierz = 2*(size(ParyObrotowe, 1) + size(ParyPostepowe, 1))+1;

%Wymuszenia postępowe
for i=1:size(WymPostepowe,1)
    FT(indeks_macierz, 1) = -1*Fab_prim(WymPostepowe(i, 1) , czas , WymPostepowe);
    indeks_macierz = indeks_macierz + 1;
end

FQ = MacierzJacobiego(Q, ParyObrotowe, ParyPostepowe, WymPostepowe);
dQ = -FQ \ FT;

end
function [d2Q] = Przyspieszenie(dQ,Q, ParyObrotowe, ParyPostepowe, WymPostepowe,czas)
%Funkcja wyznaczająca przyśpieszenie

omega = [0 -1; 1 0];
gamma = zeros(size(ParyObrotowe,1)+size(ParyPostepowe,1)+size(WymPostepowe,1), 1);
indeks_macierzy = 1;

%Pary obrotowe
for i=1:size(ParyObrotowe,1)
    if(ParyObrotowe(i, 1) == 0) %Dla podstawy
        Roti = Rot(0);
        dFii=0;
    else
        Roti = Rot(Q(3*ParyObrotowe(i,1)));
        dFii=dQ(3*ParyObrotowe(i,1));
    end
    if(ParyObrotowe(i, 2) == 0) %Dla podstawy
        Rotj = Rot(0);
        dFij=0;
    else
        Rotj = Rot(Q(3*ParyObrotowe(i,2)));
        dFij=dQ(3*ParyObrotowe(i,2));
    end

    %Wzór 2.42 ze skryptu
    gamma(indeks_macierzy:indeks_macierzy+1, 1) = Roti * ParyObrotowe(i, 3:4)' * dFii^2 - Rotj * ParyObrotowe(i, 5:6)' * dFij^2;

    indeks_macierzy = indeks_macierzy + 2;
end

%Pary postępowe
for i=1:size(ParyPostepowe,1)
    %Wzór 2.48 ze skryptu
    gamma(indeks_macierzy, 1) = 0;
    indeks_macierzy = indeks_macierzy + 1;

    if(ParyPostepowe(i, 1) == 0) %Dla podstawy
        Roti = Rot(0);
        ri=[0;0];
        dFii=0;
        dri=[0;0];
    else
        Roti = Rot(Q(3*ParyPostepowe(i,1)));
        ri=Q(3*ParyPostepowe(i,1)-2:3*ParyPostepowe(i,1)-1);
        dFii=dQ(3*ParyPostepowe(i,1));
        dri=dQ(3*ParyPostepowe(i,1)-2:3*ParyPostepowe(i,1)-1);
    end
    if(ParyPostepowe(i, 2) == 0) %Dla podstawy
        Rotj = Rot(0);
        rj=[0;0];
        dFij=0;
        drj=[0;0];
    else
        Rotj = Rot(Q(3*ParyPostepowe(i,2)));
        rj=Q(3*ParyPostepowe(i,2)-2:3*ParyPostepowe(i,2)-1);
        dFij=dQ(3*ParyPostepowe(i,2));
        drj=dQ(3*ParyPostepowe(i,2)-2:3*ParyPostepowe(i,2)-1);
    end

    %Wzór 2.59 ze skryptu
    gamma(indeks_macierzy, 1) = (Rotj * [-ParyPostepowe(i, 5); ParyPostepowe(i, 4)])'*(2 * omega * (drj - dri) * dFij + (rj - ri) * dFij^2 - Roti * ParyPostepowe(i, 6:7)' * (dFij - dFii)^2 );

    indeks_macierzy = indeks_macierzy + 1;
end

%Tu należałoby wstawić wymuszenia obrotowe

for ii=1:size(WymPostepowe,1) %petla po wszystkich wymuszeniach postepowych
    i=WymPostepowe(ii, 1);
    if(ParyPostepowe(i, 1) == 0) %Dla podstawy
        Roti = Rot(0);
        ri=[0;0];
        dFii=0;
        dri=[0;0];
    else
        Roti = Rot(Q(3*ParyPostepowe(i,1)));
        ri=Q(3*ParyPostepowe(i,1)-2:3*ParyPostepowe(i,1)-1);
        dFii=dQ(3*ParyPostepowe(i,1));
        dri=dQ(3*ParyPostepowe(i,1)-2:3*ParyPostepowe(i,1)-1);
    end
    if(ParyPostepowe(i, 2) == 0) %Dla podstawy
        Rotj = Rot(0);
        rj=[0;0];
        dFij=0;
        drj=[0;0];
    else
        Rotj = Rot(Q(3*ParyPostepowe(i,2)));
        rj=Q(3*ParyPostepowe(i,2)-2:3*ParyPostepowe(i,2)-1);
        dFij=dQ(3*ParyPostepowe(i,2));
        drj=dQ(3*ParyPostepowe(i,2)-2:3*ParyPostepowe(i,2)-1);
    end

    %Wzór 2.59 ze skryptu
    gamma(indeks_macierzy, 1) = (Rotj * ParyPostepowe(i, 4:5)')'*(2 * omega * (drj - dri) * dFij + (rj - ri) * dFij^2 - Roti * ParyPostepowe(i, 6:7)' * (dFij - dFii)^2 ) +Fab_bis(i, czas, WymPostepowe);

    indeks_macierzy = indeks_macierzy + 1;
end

FQ = MacierzJacobiego(Q, ParyObrotowe, ParyPostepowe, WymPostepowe);
d2Q = FQ \ gamma;
end
function [q] = PunktPolozenie(Q,nr_czlonu,S_A)
q=Q(3*nr_czlonu-2:3*nr_czlonu-1,:);
for i=1:length(q)
    q(:,i)=q(:,i)+Rot(Q(3*nr_czlonu,i))*S_A;
end
end
function [dq] = PunktPredkosc(Q,dQ,nr_czlonu,S_A)
omega = [0 -1; 1 0];
dq=dQ(3*nr_czlonu-2:3*nr_czlonu-1,:);
for i=1:length(dq)
    dq(:,i)=dq(:,i)+(omega*Rot(Q(3*nr_czlonu,i))*S_A*dQ(3*nr_czlonu,i));
end
end
function [dq] = PunktPrzyspieszenie(Q,dQ,d2Q,nr_czlonu,S_A)
omega = [0 -1; 1 0];
dq=d2Q(3*nr_czlonu-2:3*nr_czlonu-1,:);
for i=1:length(dq)
    dq(:,i)=dq(:,i)+(omega*Rot(Q(3*nr_czlonu,i))*S_A*d2Q(3*nr_czlonu,i))-(Rot(Q(3*nr_czlonu,i))*S_A*dQ(3*nr_czlonu,i)^2);
end
end
function R = Rot(fi)
%Macierz rotacji dla przypadku 2D - kąt fi wyrażony w radianach
R = [cos(fi) -sin(fi); sin(fi) cos(fi)];  
end
function d2q = Uklad(Czlony,dCzlony, ParyObrotowe, ParyPostepowe, WymPost, Masy,Sprezyto_tlumiace,Sily)
alfa_bombardier=10;
beta_bombardier=10;

Jakobi = MacierzJacobiego(Czlony, ParyObrotowe, ParyPostepowe, WymPost);
Q=WektorSil(Masy,Sily,Czlony,dCzlony,Sprezyto_tlumiace,ParyPostepowe);
wekt_gamma=Wektor_Gamma(dCzlony,Czlony,ParyObrotowe,ParyPostepowe,WymPost);

A=[Masy Jakobi'; Jakobi zeros(size(Jakobi,1))];
b=[Q;wekt_gamma-2*alfa_bombardier*Jakobi*dCzlony-beta_bombardier^2*Wiezy(Czlony,ParyObrotowe,ParyPostepowe,WymPost)];

x=A\b;
d2q=x(size(x,1)/2+1:size(x,1));
end
% function [ret] = Uklad_ode(t,Y, ParyObrotowe, ParyPostepowe, WymPost, Masy,Sprezyto_tlumiace,Sily)
% d2Q=Uklad(Y(1:size(Masy,1)),Y(size(Masy,1)+1:size(Y,1)), ParyObrotowe, ParyPostepowe, WymPost, Masy,Sprezyto_tlumiace,Sily);
% ret=[Y(size(Masy,1)+1:size(Y,1));d2Q];
% end
function [ret] = Uklad_ode(t,Y, ParyObrotowe, ParyPostepowe, WymPost, Masy,Sprezyto_tlumiace,Sily)
alfa_bombardier=10;
beta_bombardier=10;

Jakobi = MacierzJacobiego(Y(1:size(Masy,1)), ParyObrotowe, ParyPostepowe, WymPost);
Q=WektorSil(Masy,Sily,Y(1:size(Masy,1)),Y(size(Masy,1)+1:size(Y,1)),Sprezyto_tlumiace,ParyPostepowe);
wekt_gamma=Wektor_Gamma(Y(size(Masy,1)+1:size(Y,1)),Y(1:size(Masy,1)),ParyObrotowe,ParyPostepowe,WymPost);

A=[Masy Jakobi'; Jakobi zeros(size(Jakobi,1))];
b=[Q;wekt_gamma-2*alfa_bombardier*Jakobi*Y(size(Masy,1)+1:size(Y,1))-beta_bombardier^2*Wiezy(Y(1:size(Masy,1)),ParyObrotowe,ParyPostepowe,WymPost)];

x=A\b;
ret=[Y(size(Masy,1)+1:size(Y,1));x(1:size(Masy,1))];
end
function [Czlony, ParyObrotowe, ParyPostepowe, WymPost, Masy, Sily, Sprezyto_tlumiace] = WczytajDane()
%     Otwarcie plików
plik_czlony = fopen('Dane/Człony.txt', 'r');
plik_pary = fopen('Dane/Pary.txt', 'r');
plik_wiezy = fopen('Dane/Wymuszenia.txt', 'r');
plik_masy = fopen('Dane/Masy.txt', 'r');
plik_sily = fopen('Dane/Sily.txt', 'r');
plik_sprezyny = fopen('Dane/Sprezysto_tlumiace.txt', 'r');

%     Wczytywanie ilości członów
ilosc_czlonow = str2double(fgetl(plik_czlony));
%     Wczytanie członów do mecierzy członów (x, y, fi)
for i = 1:ilosc_czlonow
    Czlony(3*i-2:3*i,1) = str2num(fgetl(plik_czlony));
end
%     Wczytanie ilości par obrotowych
ilosc_par_obr = str2double(fgetl(plik_pary));
ParyObrotowe = zeros(ilosc_par_obr, 6);
%     Punkt połączenia pary kinematycznej
A = [0;0];
fi = [0 0];
for i = 1:ilosc_par_obr
    linijka = str2num(fgetl(plik_pary));
    [ParyObrotowe(i,1), ParyObrotowe(i,2), A(1), A(2)] = deal(linijka(1), linijka(2), linijka(3), linijka(4));

    for j = 1:2
        %             Try aby wyłapać człon utwierdzony
        try
            fi(j) = Czlony(3*ParyObrotowe(i,j));
            %                 ParyObrotowe(1,3:4) = Rot(fi1)'*(A-Czlon1(x,y))
            ParyObrotowe(i,j*2+1:j*2+2) = Rot(fi(j))'*(A-Czlony(3*ParyObrotowe(i,j)-2 : 3*ParyObrotowe(i,j)-1));
        catch
            ParyObrotowe(i,j*2+1:j*2+2) = A;
        end
    end
end
% Wczytanie ilości par postępowych
ilosc_par_post = str2double(fgetl(plik_pary));
ParyPostepowe = zeros(ilosc_par_post,10);
B = [0;0];
for i = 1:ilosc_par_post
    linijka = str2num(fgetl(plik_pary));
    [ParyPostepowe(i,1),ParyPostepowe(i,2),A(1),A(2),B(1),B(2)] = deal(linijka(1),linijka(2),linijka(3),linijka(4),linijka(5),linijka(6));
    for j = 1:2
        %         Try aby wyłapać człon utwierdzony
        try
            fi(j) = Czlony(3*ParyPostepowe(i,j));
            % ParyObrotowe(1,6:7) = Rot(fi1)'*(A-Czlon1(x,y))
            if(j==1)
                ParyPostepowe(i,j*2+4:j*2+5) = Rot(fi(j))'*(A-Czlony(3*ParyPostepowe(i,j)-2 : 3*ParyPostepowe(i,j)-1));
            else
                ParyPostepowe(i,j*2+4:j*2+5) = Rot(fi(j))'*(B-Czlony(3*ParyPostepowe(i,j)-2 : 3*ParyPostepowe(i,j)-1));
            end
        catch
            fi(j) = 0;
            if(j==1)
                ParyPostepowe(i,j*2+4:j*2+5) = A;
            else
                ParyPostepowe(i,j*2+4:j*2+5) = B;
            end
        end
    end
    ParyPostepowe(i,3) = fi(1) - fi(2);
    U = [B(1)-A(1);B(2)-A(2)];
    ParyPostepowe(i,10)=norm(U);
    U = U/norm(U);
    U = Rot(fi(2))'*U;
    ParyPostepowe(i,4:5) = U';

end

% Wymuszenia postępowe
ilosc_wym_post = str2num(fgetl(plik_wiezy));
WymPost = zeros(ilosc_wym_post,4);
for i  = 1:ilosc_wym_post
    WymPost(i,:) = str2num(fgetl(plik_wiezy));
end

% Macierz masowa
Masy=zeros(ilosc_czlonow,ilosc_czlonow);
for i = 1:ilosc_czlonow
    linijka = str2num(fgetl(plik_masy));
    Masy(3*i-2, 3*i-2) = linijka(2);
    Masy(3*i-1, 3*i-1) = linijka(2);
    Masy(3*i, 3*i) = linijka(3);
end

% Sily zewnetrzne (oprocz grawitacji)
ilosc_sil = str2num(fgetl(plik_sily));
Sily=zeros(ilosc_sil,5);
for i = 1:ilosc_sil
    linijka = str2num(fgetl(plik_sily));
    Sily(i,1)=linijka(1);
    Sily(i,2:3)=Rot(deg2rad(linijka(3)))*[linijka(2) 0]';
    Sily(i,4:5)=Rot(Czlony(Sily(i,1)*3))*(linijka(4:5)'-Czlony(Sily(i,1)*3-2:Sily(i,1)*3-1));
end

% Elementy sprezysto-tlumiace
ilosc_sprezyn = str2num(fgetl(plik_sprezyny));
Sprezyto_tlumiace=zeros(ilosc_sprezyn,3+4);
for i =1:ilosc_sprezyn
    Sprezyto_tlumiace(i,1:3+4) = str2num(fgetl(plik_sprezyny));
    nr_pary=Sprezyto_tlumiace(i,1);
    czlon_1=Czlony(ParyPostepowe(nr_pary,1)*3-[2;1]);
    czlon_2=Czlony(ParyPostepowe(nr_pary,2)*3-[2;1]);
    Sprezyto_tlumiace(i,4:5) = czlon_2-Sprezyto_tlumiace(i,4:5)'; %[-0.15;-0.45];
    Sprezyto_tlumiace(i,6:7) = czlon_1-Sprezyto_tlumiace(i,6:7)'; %[0.15;0.45];
end
end
function [gamma] = Wektor_Gamma(dQ,Q, ParyObrotowe, ParyPostepowe, WymPostepowe)
%Funkcja wyznaczająca przyśpieszenie

omega = [0 -1; 1 0];
gamma = zeros(size(ParyObrotowe,1)+size(ParyPostepowe,1)+size(WymPostepowe,1), 1);
indeks_macierzy = 1;

%Pary obrotowe
for i=1:size(ParyObrotowe,1)
    if(ParyObrotowe(i, 1) == 0) %Dla podstawy
        Roti = Rot(0);
        dFii=0;
    else
        Roti = Rot(Q(3*ParyObrotowe(i,1)));
        dFii=dQ(3*ParyObrotowe(i,1));
    end
    if(ParyObrotowe(i, 2) == 0) %Dla podstawy
        Rotj = Rot(0);
        dFij=0;
    else
        Rotj = Rot(Q(3*ParyObrotowe(i,2)));
        dFij=dQ(3*ParyObrotowe(i,2));
    end

    %Wzór 2.42 ze skryptu
    gamma(indeks_macierzy:indeks_macierzy+1, 1) = Roti * ParyObrotowe(i, 3:4)' * dFii^2 - Rotj * ParyObrotowe(i, 5:6)' * dFij^2;

    indeks_macierzy = indeks_macierzy + 2;
end

%Pary postępowe
for i=1:size(ParyPostepowe,1)
    %Wzór 2.48 ze skryptu
    gamma(indeks_macierzy, 1) = 0;
    indeks_macierzy = indeks_macierzy + 1;

    if(ParyPostepowe(i, 1) == 0) %Dla podstawy
        Roti = Rot(0);
        ri=[0;0];
        dFii=0;
        dri=[0;0];
    else
        Roti = Rot(Q(3*ParyPostepowe(i,1)));
        ri=Q(3*ParyPostepowe(i,1)-2:3*ParyPostepowe(i,1)-1);
        dFii=dQ(3*ParyPostepowe(i,1));
        dri=dQ(3*ParyPostepowe(i,1)-2:3*ParyPostepowe(i,1)-1);
    end
    if(ParyPostepowe(i, 2) == 0) %Dla podstawy
        Rotj = Rot(0);
        rj=[0;0];
        dFij=0;
        drj=[0;0];
    else
        Rotj = Rot(Q(3*ParyPostepowe(i,2)));
        rj=Q(3*ParyPostepowe(i,2)-2:3*ParyPostepowe(i,2)-1);
        dFij=dQ(3*ParyPostepowe(i,2));
        drj=dQ(3*ParyPostepowe(i,2)-2:3*ParyPostepowe(i,2)-1);
    end

    %Wzór 2.59 ze skryptu
    gamma(indeks_macierzy, 1) = (Rotj * [-ParyPostepowe(i, 5); ParyPostepowe(i, 4)])'*(2 * omega * (drj - dri) * dFij + (rj - ri) * dFij^2 - Roti * ParyPostepowe(i, 6:7)' * (dFij - dFii)^2 );

    indeks_macierzy = indeks_macierzy + 1;
end
end
function [Q]=WektorSil(Masy,Sily,Czlony,dCzlony,Sprezyto_tlumiace,ParyPostepowe)
g=[0 -9.81]';
omega = [0 -1; 1 0];

% grawitacja
Q_g=zeros(length(Masy),1);
for i=1:length(Masy)/3
    Q_g(3*i-2:3*i-1)=Masy(3*i-2,3*i-2)*g;
end

% sily zewnetrzne
Q_fz=zeros(length(Masy),1);
for i=1:size(Sily)
    id_czlonu=Sily(i,1);
    Q_fz(id_czlonu*3-2:id_czlonu*3-1)=Sily(i,2:3);
    Q_fz(id_czlonu*3)=(omega*Rot(Czlony(3*id_czlonu))*Sily(i,4:5)')'*Sily(i,2:3)';
end

% sily sprezystosci
Q_spr=zeros(length(Masy),1);
for iter=1:1:size(Sprezyto_tlumiace,1)
    c2=Czlony(ParyPostepowe(iter,2)*3-[2 1]);
    c1=Czlony(ParyPostepowe(iter,1)*3-[2 1]);
    v_c1=dCzlony(ParyPostepowe(iter,1)*3-[2 1]);
    v_c2=dCzlony(ParyPostepowe(iter,2)*3-[2 1]);
    omega_c1=dCzlony(ParyPostepowe(iter,1)*3);
    omega_c2=dCzlony(ParyPostepowe(iter,2)*3);
    fi1=Czlony(ParyPostepowe(iter,1)*3);
    fi2=Czlony(ParyPostepowe(iter,2)*3);
    d = c2 - c1 + Rot(fi2)*Sprezyto_tlumiace(iter,6:7)' - Rot(fi1)*Sprezyto_tlumiace(iter,4:5)';
    u = d/norm(d);
    
    Fk = Sprezyto_tlumiace(iter,2)*(norm(d) - ParyPostepowe(iter,10));
    Fc = Sprezyto_tlumiace(iter,3) * u' * (v_c2 - v_c1 + ...
        omega * Rot(fi2) * Sprezyto_tlumiace(iter,6:7)' * omega_c2 - ...
        omega * Rot(fi1) * Sprezyto_tlumiace(iter,4:5)' * omega_c1);
    
    if(ParyPostepowe(iter,1)>0)
        Q_spr(3*ParyPostepowe(iter,1)+(-2:0)) = Q_spr(3*ParyPostepowe(iter,1)+(-2:0)) + ...
            [eye(2); (omega * Rot(fi1)*Sprezyto_tlumiace(iter,4:5)')'] * u * (Fc+Fk);
    end
    if(ParyPostepowe(iter,2)>0)
        Q_spr(3*ParyPostepowe(iter,2)+(-2:0)) = Q_spr(3*ParyPostepowe(iter,2)+(-2:0)) + ...
            [eye(2); (omega * Rot(fi2)*Sprezyto_tlumiace(iter,6:7)')'] * (-u) * (Fc+Fk);
    end
end

Q=Q_spr+Q_fz+Q_g;
end
function F = Wiezy(Q, ParyObrotowe, ParyPostepowe, WymPostepowe)
F = zeros(size(ParyObrotowe,1)+size(ParyPostepowe,1)+size(WymPostepowe,1),1);
indeks_macierz = 1;
% Pary obrotowe
for i = 1:(size(ParyObrotowe,1))
    if ParyObrotowe(i,1)~=0
        ri=Q(3*ParyObrotowe(i,1)-2:3*ParyObrotowe(i,1)-1);
        Roti=Rot(Q(3*ParyObrotowe(i,1)));
    else
        ri=[0;0];
        Roti=Rot(0);
    end
    if ParyObrotowe(i,2)~=0
        rj=Q(3*ParyObrotowe(i,2)-2:3*ParyObrotowe(i,2)-1);
        Rotj=Rot(Q(3*ParyObrotowe(i,2)));
    else
        rj=[0;0];
        Rotj=Rot(0);
    end
    % Wzór 2.18 w skrypcie z wykładu
    F(indeks_macierz:indeks_macierz+1,1) = ri + Roti * ParyObrotowe(i,3:4)' - (rj + Rotj * ParyObrotowe(i,5:6)');
    indeks_macierz = indeks_macierz + 2;
end

% Pary postepowe
for i = 1:(size(ParyPostepowe,1))
    if ParyPostepowe(i,1)~=0
        ri=Q(3*ParyPostepowe(i,1)-2:3*ParyPostepowe(i,1)-1);
        Roti=Rot(Q(3*ParyPostepowe(i,1)));
        fii=Q(3*ParyPostepowe(i,1));
    else
        ri=[0;0];
        Roti=Rot(0);
        fii=0;
    end
    if ParyPostepowe(i,2)~=0
        rj=Q(3*ParyPostepowe(i,2)-2:3*ParyPostepowe(i,2)-1);
        Rotj=Rot(Q(3*ParyPostepowe(i,2)));
        fij=Q(3*ParyPostepowe(i,2));
    else
        rj=[0;0];
        Rotj=Rot(0);
        fij=0;
    end
    %Wzór 2.19
    F(indeks_macierz,1) = fii - fij - ParyPostepowe(i,3);
    %Wzór 2.22
    F(indeks_macierz+1,1) = (Rotj * [-ParyPostepowe(i,5); ParyPostepowe(i,4)])'*(rj - ri - Roti * ParyPostepowe(i,6:7)') + [-ParyPostepowe(i,5); ParyPostepowe(i,4)]' * ParyPostepowe(i,8:9)';

    indeks_macierz = indeks_macierz + 2;
end

% % wymuszenia postepowe
% for ii = 1:(size(WymPostepowe,1))
%     i=WymPostepowe(ii,1);
%     if ParyPostepowe(i,1)~=0
%         ri=Q(3*ParyPostepowe(i,1)-2:3*ParyPostepowe(i,1)-1);
%         Roti=Rot(Q(3*ParyPostepowe(i,1)));
%     else
%         ri=[0;0];
%         Roti=Rot(0);
%     end
%     if ParyPostepowe(i,2)~=0
%         rj=Q(3*ParyPostepowe(i,2)-2:3*ParyPostepowe(i,2)-1);
%         Rotj=Rot(Q(3*ParyPostepowe(i,2)));
%     else
%         rj=[0;0];
%         Rotj=Rot(0);
%     end
% 
%     % wzor 2.28 na wymuszenie w parze postepowej
%     F(indeks_macierz,1) = (Rotj * ParyPostepowe(i,4:5)')'*(rj + Rotj * ParyPostepowe(i,8:9)' - ri - Roti * ParyPostepowe(i,6:7)') - Fab(i,czas,ParyPostepowe, WymPostepowe);
% 
%     indeks_macierz = indeks_macierz + 1;
% end

end
function x = Fab(id_wiezu, t,ParyPostepowe, WymPostepowe)
% wektorowa funkcja czasu
%lk=sqrt((ParyPostepowe(id_wiezu,10)-ParyPostepowe(id_wiezu,12))^2+(ParyPostepowe(id_wiezu,11)-ParyPostepowe(id_wiezu,13))^2);
A=WymPostepowe(id_wiezu,2);
omega=WymPostepowe(id_wiezu,3);
fi=WymPostepowe(id_wiezu,4);

x = ParyPostepowe(id_wiezu,10) + A*sin(omega*t + fi);
end
function x = Fab_bis(id_wiezu, t, WymPostepowe)
% wektorowa funkcja czasu
A=WymPostepowe(id_wiezu,2);
omega=WymPostepowe(id_wiezu,3);
fi=WymPostepowe(id_wiezu,4);

x = -A*omega*omega*sin(omega*t + fi);,
end
function x = Fab_prim(id_wiezu, t, WymPostepowe)
% wektorowa funkcja czasu
A=WymPostepowe(id_wiezu,2);
omega=WymPostepowe(id_wiezu,3);
fi=WymPostepowe(id_wiezu,4);

x = A*omega*cos(omega*t + fi);
end
