% ======= Start.m (kod wykonywalny) =======
DT = 0.1
T_MAX = 30
[T,Q, DQ, DDQ] = Main(DT,T_MAX);
plot(T, Q(2, :))
xlabel('Czas [s]')
ylabel('Pozycja Y członu c1 [m]')
title('Pozycja Y członu c1 w czasie')
grid on
%plot(T,Q(n,:))

% czlon a
% a(x) - n=3*a-2
% a(y) - n=3*a-1
% a(fi) - n=3*a


% ======= Main.m =======
function [T,Q, DQ, DDQ] = Main(DT,T_MAX)
    
    %Pobranie danych z plików txt
    [Czlony, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe] = Dane();
    
    % Przybliżenie startowe  
    Q_TEMP = Czlony;
    DQ_TEMP = zeros(length(Q_TEMP),1); 
    DDQ_TEMP = zeros(length(Q_TEMP),1);  
    
    %Wyznaczenie liczby wyników i zadeklarowanie licznika 
    IleIteracji = 1;
    Size = T_MAX/DT + 1;
    
    %Zadeklarowanie rozmiaru macierzy i wektorów gromadzących wyniki
    T = zeros(1, Size);
    Q = zeros(length(Q_TEMP), Size);
    DQ = zeros(length(Q_TEMP), Size);
    DDQ = zeros(length(Q_TEMP), Size);

    %Rozwiązywanie zadań kinematyki w kolejnych chwilach T_TEMP
    for T_TEMP = 0:DT:T_MAX
        Q_TEMP=Q_TEMP+DQ_TEMP*DT+0.5*DDQ_TEMP*DT*DT;
        
        %Rozwiązanie zadania o położenie
        Q_TEMP=NewRaph(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe, T_TEMP); 
        
        %Rozwiązanie zadania o prędkości
        DQ_TEMP=Predkosc(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe,T_TEMP); 
        
        %Rozwiązanie zadania o przyśpieszenia 
        DDQ_TEMP=Przyspieszenie(DQ_TEMP,Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe,T_TEMP);

        %Zapisane wyników do macierzy i wektorów gromadzących
        T(1,IleIteracji) = T_TEMP; 
        Q(:,IleIteracji) = Q_TEMP;
        DQ(:,IleIteracji) = DQ_TEMP;
        DDQ(:,IleIteracji) = DDQ_TEMP;
    
        IleIteracji = IleIteracji + 1;
    end
end


% ======= Dane.m =======
function [Czlony, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe] = Dane()
    %Funkcja wczytująca dane z plików txt 
    
    c = fopen('Dane/Człony.txt', 'r');
    p = fopen('Dane/Pary.txt', 'r');
    w = fopen('Dane/Wymuszenia.txt', 'r');
    
    %Wczytanie członów
    no = str2num(fgetl(c));  % wczytanie liczby członów
    for i = 1:no
        Czlony(3*(i-1)+1:3*i) = str2num(fgetl(c));  %wczytanie poszczególnych członów [x,y,fi]
    end
    Czlony
    Czlony = Czlony'   %macierz 3 x no
    
    %Wczytanie par obrotowych - rozdzial 2.3
    no = str2num(fgetl(p));     %pierwsza dana to liczba par obrotowych
    ParyObrotowe = zeros(no, 6);    
    for i = 1:no
        tmp = str2num(fgetl(p));
        Czlon1 = tmp(1); Czlon2 = tmp(2);
        Polozenie = [tmp(3); tmp(4)];

        if(Czlon1==0)
            q1 = [0;0]; f1 = 0;
        else
            q1 = Czlony(3*Czlon1-2 : 3*Czlon1-1); f1 = Czlony(3*Czlon1);
        end
        
        if(Czlon2==0)
            q2 = [0;0]; f2 = 0;
        else
            q2 = Czlony(3*Czlon2-2 : 3*Czlon2-1); f2 = Czlony(3*Czlon2);
        end
        
        R1 = Rot(f1); R2 = Rot(f2);

        ParyObrotowe(i,1:2) = [Czlon1,Czlon2];
        ParyObrotowe(i,3:4) = R1'*(Polozenie-q1);
        ParyObrotowe(i,5:6) = R2'*(Polozenie-q2);
        
    end
    ParyObrotowe
    %Wczytanie par postępowych - rozdział 2.3
    no = str2num(fgetl(p));
    ParyPostepowe = zeros(no,13);
    for i = 1:no
        tmp = str2num(fgetl(p));
        Czlon1 = tmp(1); Czlon2 = tmp(2);
        Polozenie1 = [tmp(3); tmp(4)];
        Polozenie2 = [tmp(5); tmp(6)];

        if(Czlon1==0)
            q1 = [0;0]; f1 = 0;
        else
            q1 = Czlony(3*Czlon1-2 : 3*Czlon1-1); f1 = Czlony(3*Czlon1);
        end
        
        if(Czlon2==0)
            q2 = [0;0]; f2 = 0;
        else
            q2 = Czlony(3*Czlon2-2 : 3*Czlon2-1); f2 = Czlony(3*Czlon2);
        end

    R1 = Rot(f1); R2 = Rot(f2);
    ParyPostepowe(i,1:2) = [Czlon1,Czlon2];
    ParyPostepowe(i,3) = f1 - f2;
    ParyPostepowe(i,6:7) = R1'*(Polozenie1-q1);
    ParyPostepowe(i,8:9) = R2'*(Polozenie2-q2);
    ParyPostepowe(i,10:13)=[Polozenie1',Polozenie2']
    U = [Polozenie2(1)-Polozenie1(1);Polozenie2(2)-Polozenie1(2)];
    U = U/norm(U);
    U = R2'*U;
    ParyPostepowe(i,4:5) = U';
    end
    ParyPostepowe
    %Wczytanie wymuszen obrotowych 
    no = str2num(fgetl(w));
    WymObrotowe = zeros(no,2);
    for i  = 1:no
        WymObrotowe(i,:) = str2num(fgetl(w));
    end
    WymObrotowe
    %Wczytanie wymuszen postepowych 
    no = str2num(fgetl(w));


    WymPostepowe = zeros(no,5);
    for i  = 1:no
        WymPostepowe(i,:) = str2num(fgetl(w));
        
    end
    WymPostepowe
end
        



% ======= NewRaph.m =======
function Q_TEMP = NewRaph(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe, T_TEMP)
    %Rozwiązanie układu równań metodą Newtona-Raphsona

    F = ones(length(Q_TEMP),1);
    EPS = 10e-9;
    IleIteracji = 1;
    MaxIteracji = 25;

    while(norm(F) > EPS && IleIteracji < MaxIteracji)
        F = Wiezy(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe, T_TEMP);
        FQ = Jakobian(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe);
        Q_TEMP = Q_TEMP - FQ\F;
        IleIteracji = IleIteracji + 1;
    end
    
    
    if(IleIteracji > MaxIteracji)
        disp('Nie udało się uzyskać zbieżności po założonej liczbie iteracji')
    end
end



% ======= Predkosc.m =======
function [DQ] = Predkosc(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe,T_TEMP)
    %Funkcja rozwiązująca zadanie o prędkości
    
    FT = zeros(length(Q_TEMP), 1);
    Pozycja = 2*(size(ParyObrotowe, 1) + size(ParyPostepowe, 1))+1;

    %Tutaj wstawiłbym wymuszenia obrotowe, ale ich nie ma
    
    %Wymuszenia postępowe
    for i=1:size(WymPostepowe,1)
        FT(Pozycja, 1) = -1*WymuszenieP(WymPostepowe(i, 2), 1, T_TEMP,ParyPostepowe, WymPostepowe);
        Pozycja = Pozycja + 1;
    end

    FQ = Jakobian(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe);
    DQ = -FQ \ FT;

    end



% ======= Przyspieszenie.m =======
function [DDQ] = Przyspieszenie(DQ_TEMP,Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe,T_TEMP)
    %Funkcja wyznaczająca przyśpieszenie

    Om = [0 -1; 1 0]; 
    Gamma = zeros(length(Q_TEMP), 1);
    Pozycja = 1;
     
    %Pary obrotowe
    for i=1:size(ParyObrotowe,1) 
        P1 = ParyObrotowe(i, 1);
        P2 = ParyObrotowe(i, 2);
    
        [~,~,Rot1] = Polozenia(Q_TEMP, P1);
        [~,~,Rot2] = Polozenia(Q_TEMP, P2);
    
        sA = ParyObrotowe(i, 3:4)';
        sB = ParyObrotowe(i, 5:6)';
    
        [~,Dfi1,~]= Polozenia(DQ_TEMP,P1);
        [~,DFi2,~]= Polozenia(DQ_TEMP,P2);
        
        %Wzór 2.40 ze skryptu
        Gamma(Pozycja:Pozycja+1, 1) = Rot1 * sA * Dfi1^2 - Rot2 * sB * DFi2^2;
    
        Pozycja = Pozycja + 2;
    end
    
    %Pary postępowe
    for i=1:size(ParyPostepowe,1) 
        %Wzór 2.46 ze skryptu        
        Gamma(Pozycja, 1) = 0;
        Pozycja = Pozycja + 1;

        P1 = ParyPostepowe(i, 1);
        P2 = ParyPostepowe(i, 2);

        U = ParyPostepowe(i, 4:5)';
        V = [-U(2); U(1)];
        sA = ParyPostepowe(i, 6:7)';

        [R1,~,Rot1] = Polozenia(Q_TEMP, P1);
        [R2,~,Rot2] = Polozenia(Q_TEMP, P2);

        [DR1,Dfi1,~] = Polozenia(DQ_TEMP,P1);
        [DR2,DFi2,~] = Polozenia(DQ_TEMP,P2);

        %Wzór 2.57 ze skryptu
        Gamma(Pozycja, 1) = (Rot2 * V)'*(2 * Om * (DR2 - DR1) * DFi2 + (R2 - R1) * DFi2^2 - Rot1 * sA * (DFi2 - Dfi1)^2 ); 

        Pozycja = Pozycja + 1;
    end

    %Tu należałoby wstawić wymuszenia obrotowe

    for i=1:size(WymPostepowe,1) %petla po wszystkich wymuszeniach postepowych    
        P1 = ParyPostepowe(WymPostepowe(i, 1), 1);
        P2 = ParyPostepowe(WymPostepowe(i, 1), 2);

        U = ParyPostepowe(WymPostepowe(i, 1), 4:5)';
        sA = ParyPostepowe(WymPostepowe(i, 1), 6:7)';

        [R1,~,Rot1] = Polozenia(Q_TEMP, P1);
        [R2,~,Rot2] = Polozenia(Q_TEMP, P2);

        [DR1,Dfi1,~]=Polozenia(DQ_TEMP,P1);
        [DR2,DFi2,~]=Polozenia(DQ_TEMP,P2);

        %Wzór 2.57 ze skryptu
        Gamma(Pozycja, 1) = (Rot2 * U)'*(2 * Om * (DR2 - DR1) * DFi2 + (R2 - R1) * DFi2^2 - Rot1 * sA * (DFi2 - Dfi1)^2 ) +WymuszenieP(WymPostepowe(i, 2), 2, T_TEMP,ParyPostepowe, WymPostepowe); 

        Pozycja = Pozycja + 1;
    end

    %Obliczenie macierzy układu równań
    FQ = Jakobian(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe);

    %Rozwiązanie
    DDQ = FQ \ Gamma;
end



% ======= Wiezy.m =======
function F = Wiezy(Q_TEMP, ParyObrotowe, ParyPostepowe, ~, WymPostepowe, T_TEMP)
    %Funkcja wyznaczająca wartości funkcji opisujących więzy
 
    %Zadeklarowanie wektora wyjściowego
    F = zeros(length(Q_TEMP),1);
    Pozycja = 1;

    %Stworzenie równań dla par obrotowych
    for i = 1:(size(ParyObrotowe,1))
        [R1, ~, Rot1] = Polozenia(Q_TEMP,ParyObrotowe(i,1));
        [R2, ~, Rot2] = Polozenia(Q_TEMP,ParyObrotowe(i,2));
        sA = ParyObrotowe(i,3:4)';
        sB = ParyObrotowe(i,5:6)';
 
        %Wzór 2.16 w skrypcie z wykładu
        F(Pozycja:Pozycja+1,1) = R1 + Rot1 * sA - (R2 + Rot2 * sB);
        
        Pozycja = Pozycja + 2;
    end

    %Stworzenie równań dla par postępowych 
    for i = 1:(size(ParyPostepowe,1))
        [R1, Fi1, Rot1] = Polozenia(Q_TEMP,ParyPostepowe(i,1));
        [R2, Fi2, Rot2] = Polozenia(Q_TEMP,ParyPostepowe(i,2));  
       
        Fi0 = ParyPostepowe(i,3);
        sA = ParyPostepowe(i,6:7)';
        sB = ParyPostepowe(i,8:9)';
        U = ParyPostepowe(i,4:5)';
        V = [-U(2); U(1)];

        %Wzór 2.17
        F(Pozycja,1) = Fi1 - Fi2 - Fi0;
        
        %Wzór 2.20
        F(Pozycja+1,1) = (Rot2 * V)'*(R2 - R1 - Rot1 * sA) + V' * sB;
        
        Pozycja = Pozycja + 2;
    end

    %Gdyby w zadaniu były wymuszenia obrotowe to powinien znaleźć się tu kod który by
    %coś z nimi robił, ale ich nie ma

    %Stworzenie równań dla wymuszeń postępowych
    for i = 1:(size(WymPostepowe,1))
        [R1,~,Rot1] = Polozenia(Q_TEMP,ParyPostepowe(WymPostepowe(i,1),1));
        [R2,~,Rot2] = Polozenia(Q_TEMP,ParyPostepowe(WymPostepowe(i,1),2));
        U = ParyPostepowe(WymPostepowe(i,1),4:5)';
        sA = ParyPostepowe(WymPostepowe(i,1),6:7)';
        sB = ParyPostepowe(WymPostepowe(i,1),8:9)';
        
        % wzor 2.26 na wymuszenie w parze postepowej
        F(Pozycja,1) = (Rot2 * U)'*(R2 + Rot2 * sB - R1 - Rot1 * sA) - WymuszenieP(WymPostepowe(i,2),0,T_TEMP,ParyPostepowe, WymPostepowe);
        
        Pozycja = Pozycja + 1;
    end
end
%% jakobian

function [FQ] = Jakobian(Q_TEMP, ParyObrotowe, ParyPostepowe, ~, WymPostepowe)
    %Funkcja wyznaczaj¹ca Jakobian równañ wiêzów 
    
    Om = [0 -1; 1 0];
    FQ = zeros(length(Q_TEMP),length(Q_TEMP));
    Pozycja = 1;

    %Pary obrotowe
    for i = 1:size(ParyObrotowe,1)
        sA = ParyObrotowe(i,3:4)';
        sB = ParyObrotowe(i,5:6)';
        [~,~,Rot1] = Polozenia(Q_TEMP,ParyObrotowe(i,1));
        [~,~,Rot2] = Polozenia(Q_TEMP,ParyObrotowe(i,2));
        
        %Je¿eli dany cz³on nie jest podstaw¹
        if ParyObrotowe(i,1) ~= 0 
            %Wzór 2.29 ze skryptu
            FQ(Pozycja:Pozycja+1, ParyObrotowe(i,1)*3-2:ParyObrotowe(i,1)*3-1) = eye(2);
            
            %Wzór 2.30 ze skryptu
            FQ(Pozycja:Pozycja+1,ParyObrotowe(i,1)*3) = Om * Rot1 * sA; 
        end
        
        %Je¿eli dany cz³on nie jest podstaw¹
        if ParyObrotowe(i,2) ~= 0
            %Wzór 2.31 ze skryptu
            FQ(Pozycja:Pozycja+1, ParyObrotowe(i,2)*3-2:ParyObrotowe(i,2)*3-1) = -eye(2); 
            
            %Wzór 2.32 ze skryptu
            FQ(Pozycja:Pozycja+1, ParyObrotowe(i,2)*3) = -Om * Rot2 * sB; 
        end
    
        Pozycja = Pozycja + 2;
    end

    %Pary postepowe
    for i = 1:size(ParyPostepowe,1)
        U = ParyPostepowe(i,4:5)';
        V = [-U(2); U(1)];
        sA = ParyPostepowe(i,6:7)';
        
        [R1,~,Rot1] = Polozenia(Q_TEMP,ParyPostepowe(i,1));
        [R2,~,Rot2] = Polozenia(Q_TEMP,ParyPostepowe(i,2));
        
        %Je¿eli dany cz³on nie jest podstaw¹
        if ParyPostepowe(i,1) ~= 0
            %Wzór 2.42 ze skryptu
            FQ(Pozycja, ParyPostepowe(i,1)*3) = 1; 
            
            %Wzór 2.47 ze skryptu
            FQ(Pozycja+1, ParyPostepowe(i,1)*3-2:ParyPostepowe(i,1)*3-1) = -(Rot2 * V)'; 
            
            %Wzór 2.48 ze skryptu
            FQ(Pozycja+1, ParyPostepowe(i,1)*3) = -(Rot2 * V)' * Om * Rot1 * sA; 
        end
        
        %Je¿eli cz³on nie jest podstaw¹
        if ParyPostepowe(i,2) ~= 0
            %Wzór 2.44 ze skryptu
            FQ(Pozycja, ParyPostepowe(i,2)*3) = -1; 
            
            %Wzór 2.49 ze skryptu
            FQ(Pozycja+1, ParyPostepowe(i,2)*3-2:ParyPostepowe(i,2)*3-1) = (Rot2 * V)'; 
            
            %Wzór 2.50 ze skryptu
            FQ(Pozycja+1, ParyPostepowe(i,2)*3) = -(Rot2 * V)' * Om * (R2 - R1 - Rot1 * sA); 
        end
   
        Pozycja = Pozycja + 2;
    end

    %Gdyby byly wymuszenia obrotowe to tutaj trzeba by by³o coœ z nimi
    %zrobiæ, ale ich nie ma

    %Wymuszenia postepowe
    for i = 1:size(WymPostepowe,1)
        U = ParyPostepowe(WymPostepowe(i,1),4:5)';
        sA = ParyPostepowe(WymPostepowe(i,1),6:7)';
        
        [R1,~,Rot1] = Polozenia(Q_TEMP,ParyPostepowe(WymPostepowe(i,1),1));
        [R2,~,Rot2] = Polozenia(Q_TEMP,ParyPostepowe(WymPostepowe(i,1),2));
        
        %Je¿eli cz³on nie jest podstaw¹
        if ParyPostepowe(WymPostepowe(i,1),1) ~= 0
            %Wzór 2.47 ze skryptu
            FQ(Pozycja, ParyPostepowe(WymPostepowe(i,1),1)*3-2:ParyPostepowe(WymPostepowe(i,1),1)*3-1) = -(Rot2 * U)'; 
            
            %Wzór 2.48 ze skryptu
            FQ(Pozycja, ParyPostepowe(WymPostepowe(i,1),1)*3) = -(Rot2 * U)' * Om * Rot1 * sA; % 
        end
        
        %Je¿eli cz³on nie jest podstaw¹
        if ParyPostepowe(WymPostepowe(i,1),2) ~= 0
            %Wzór 2.49 ze skryptu
            FQ(Pozycja, ParyPostepowe(WymPostepowe(i,1),2)*3-2:ParyPostepowe(WymPostepowe(i,1),2)*3-1) = (Rot2 * U)'; 
            
            %Wzór 2.50 ze skryptu
            FQ(Pozycja, ParyPostepowe(WymPostepowe(i,1),2)*3) = -(Rot2 * U)' * Om * (R2 - R1 - Rot1 * sA); 
        end
        Pozycja = Pozycja + 1;
    end
    
    %Sprawdzenie osobliwoœci Jakobianu
    if(det(FQ)==0)
        disp('Uwaga - Jakobian jest osobliwy :c')
    end
end
function [ x ] = WymuszenieP( indeks, pochodna, t,ParyPostepowe, WymPostepowe )
    %Zwraca dla danego czasu wartoœæ danej pochodnej wybranego wymuszenia
    %postêpowego
    
    
     L=sqrt((ParyPostepowe(indeks,10)-ParyPostepowe(indeks,12))^2+(ParyPostepowe(indeks,11)-ParyPostepowe(indeks,13))^2);
     A=WymPostepowe(indeks,3);
     OM=WymPostepowe(indeks,4);
     FI=WymPostepowe(indeks,5);

    if(pochodna == 0)
    
            x = L + A*sin(OM*t + FI);
            
       
    elseif(pochodna == 1)
        
             x = A*OM*cos(OM*t + FI);
               
    elseif(pochodna == 2)
             
             x = -A*OM*OM*sin(OM*t + FI);
                
 
    else
        disp('Zbyt du¿y b¹dŸ nieprawid³owy indeks pochodnej')
    end


end
% 
%% rotacja
function R = Rot(fi)
%Macierz rotacji dla przypadku 2D - k¹t fi wyra¿ony w radianach
R = [cos(fi) -sin(fi); sin(fi) cos(fi)];  
end
%% polozenia

function [R, Fi, Roti] = Polozenia(Q,indeks)
    %Funkcja pobieraj¹ca aktualne wspó³rzêdne 

    if(indeks == 0) %Dla podstawy
        R = [0;0]; Fi = 0; Roti = Rot(Fi);
    else
        R = Q(3*indeks-2:3*indeks-1);
        Fi = Q(3*indeks);
        Roti = Rot(Fi);
    end
end



