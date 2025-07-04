%START.m
%fprintf('Krok czasowy:\n')
dT = 0.1; %input("Wpisz dt:\n");
%fprintf('Czas symulacji [s]:\n')
T_Max = 30; %input("Wpisz t maksymalne:\n")

[T,Q, DQ, D2Q] = Main(dT,T_Max);
[Czlony,Punkty,~,~,~,~]=Dane();
[pq,pdq,pd2q] = PunktPolozenie(Punkty,Czlony,Q,DQ,D2Q,10,9); %(~,~,~,~,~,nr_czlonu,nr_punktu)
%[czlon,punkt]
%[1,1]-punkt D
%[1,2]-punkt F
%[1,3]-punkt G
%[2,4]-punkt H
%[8,5]-punkt J
%[8,6]-punkt J
%[8,7]-punkt K
%[10,8]-punkt L
%[10,9]-punkt M
%[10,10]-punkt D

%pq,pdq,pd2q - 1-x,2-y
figure()
plot(T,pq(1,:))
xlabel('Czas symulacji [s]')
ylabel('Położenie x punktu M [m]')
grid on

function R = Rot(fi)
%Funkcja obliczająca macierz rotacji2D
R = [cos(fi) -sin(fi); sin(fi) cos(fi)];  
end

%Main.m
function [T,Q, DQ, D2Q] = Main(dT,T_Max)   
    %Wczytanie parametrów z plików
    [Czlony, ~, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe] = Dane();
    
    %Wyznaczenie liczby wyników i zadeklarowanie licznika 
    IleIteracji = 1;
    Size = T_Max/dT + 1;

    % Przybliżenie startowe  
    Q_TEMP = Czlony;
    DQ_TEMP = zeros(length(Q_TEMP),1);
    D2Q_TEMP = zeros(length(Q_TEMP),1);
    
    %Zadeklarowanie rozmiaru macierzy i wektorów gromadzących wyniki
    T = zeros(1, Size);
    Q = zeros(length(Q_TEMP), Size);
    DQ = zeros(length(Q_TEMP), Size);
    D2Q = zeros(length(Q_TEMP), Size);

    %Rozwiązywanie zadań kinematyki w kolejnych chwilach T_TEMP
    for T_TEMP = 0:dT:T_Max
        Q_TEMP=Q_TEMP+DQ_TEMP*dT+0.5*D2Q_TEMP*dT*dT;
        
        %Rozwiązanie zadania o położenie
        Q_TEMP=NewRaph(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe, T_TEMP);
        %Rozwiązanie zadania o prędkości
        DQ_TEMP=Predkosc(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe,T_TEMP);
        %Rozwiązanie zadania o przyśpieszenia 
        D2Q_TEMP=Przyspieszenie(DQ_TEMP,Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe,T_TEMP);
        %Zapisane wyników do macierzy i wektorów gromadzących
        T(1,IleIteracji) = T_TEMP; 
        Q(:,IleIteracji) = Q_TEMP;
        DQ(:,IleIteracji) = DQ_TEMP;
        D2Q(:,IleIteracji) = D2Q_TEMP;
        IleIteracji = IleIteracji + 1;
    end
    %%
end
function [ x ] = WymuszenieP( indeks, pochodna, t,ParyPostepowe, WymPostepowe )
    %Zwraca dla danego czasu wartość danej pochodnej wybranego wymuszenia postępowego
    
    
     L=sqrt((ParyPostepowe(indeks,10)-ParyPostepowe(indeks,12))^2+(ParyPostepowe(indeks,11)-ParyPostepowe(indeks,13))^2)+WymPostepowe(indeks,3); %początkowa długość l_k wynika z odległości pomiędzy punktami tworzącymi kinematyczną parę postępową. Do wartości tej dodawana jest zczytywana z pliku wartość wychylenia początkowego
     A=WymPostepowe(indeks,4);
     OM=WymPostepowe(indeks,5);
     FI=WymPostepowe(indeks,6);

    if(pochodna == 0)
    
            x = L + A*sin(OM*t + FI);
            
       
    elseif(pochodna == 1)
        
             x = A*OM*cos(OM*t + FI);
               
    elseif(pochodna == 2)
             
             x = -A*OM*OM*sin(OM*t + FI);
                
 
    else
        disp('Nieprawidłowy indeks pochodnej')
    end


end


%Dane.m
function [Czlony, Punkty, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe,Masy,Sily,Sprezyto_tlumiace] = Dane()
%%
    %Funkcja wczytująca parametry z plików 
    c = fopen('Dane/Człony.txt', 'r');
    p = fopen('Dane/Pary.txt', 'r');
    w = fopen('Dane/Wymuszenia.txt', 'r');
    po=fopen('Dane/Punkty.txt','r');
    
    %Wczytanie ilości członów
    no_cz = str2num(fgetl(c));  % wczytanie liczby członów
    for i = 1:no_cz
        Czlony(3*(i-1)+1:3*i) = str2num(fgetl(c));  %wczytanie poszczególnych członów do macierzy członów [x,y,fi]
    end
    %fprintf('Czlony:\n')
    %disp(Czlony)
    Czlony = Czlony'; %macierz 3 x no

    %wczytanie punktów
     no_po = str2num(fgetl(po));  % wczytanie liczby członów
    for i = 1:no_po
        Punkty(4*(i-1)+1:4*i) = str2num(fgetl(po));  %wczytanie poszczególnych członów do macierzy członów [x,y,fi]
    end
    %fprintf('Punkty:\n')
    %disp(Punkty)
    Punkty = Punkty'; %macierz 3 x no
    %%
    %Wczytanie par obrotowych
    no_rev = str2num(fgetl(p));     %liczba par obrotowych
    ParyObrotowe = zeros(no_rev, 6); %stworzenie pustej początkowej macierzy przechowywującej pary obrotowe
    for i = 1:no_rev %przepisanie wartości z pliku wsadowego do kompilatora
        tmp = str2num(fgetl(p));
        Czlon1 = tmp(1); Czlon2 = tmp(2);
        Polozenie = [tmp(3); tmp(4)];

        if(Czlon1==0) %sprawdzenie warunku podłoża
            q1 = [0;0]; f1 = 0;
        else
            q1 = Czlony(3*Czlon1-2 : 3*Czlon1-1); f1 = Czlony(3*Czlon1);
        end
        
        if(Czlon2==0)%sprawdzenie warunku podłoża
            q2 = [0;0]; f2 = 0;
        else
            q2 = Czlony(3*Czlon2-2 : 3*Czlon2-1); f2 = Czlony(3*Czlon2);
        end
        
        R1 = Rot(f1); 
        R2 = Rot(f2);

        ParyObrotowe(i,1:2) = [Czlon1,Czlon2]; %definicja wektora par obrotowych
        ParyObrotowe(i,3:4) = R1'*(Polozenie-q1);
        ParyObrotowe(i,5:6) = R2'*(Polozenie-q2);
        
    end
    %fprintf('Pary obrotowe\n')
    %disp(ParyObrotowe)
    %Wczytanie par postępowych - analogiczne do par obrotowych
    no_pp = str2num(fgetl(p));
    ParyPostepowe = zeros(no_pp,13);
    for i = 1:no_pp
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
    ParyPostepowe(i,10:13)=[Polozenie1',Polozenie2'];
    U = [Polozenie2(1)-Polozenie1(1);Polozenie2(2)-Polozenie1(2)];
    U = U/norm(U);
    U = R2'*U;
    ParyPostepowe(i,4:5) = U';
    end
    %fprintf('Pary Postępowe:\n')
    %disp(ParyPostepowe)
    %Wczytanie wymuszen obrotowych 
    no_wymo = str2num(fgetl(w));
    WymObrotowe = zeros(no_wymo,2);
    for i  = 1:no_wymo
        WymObrotowe(i,:) = str2num(fgetl(w));
    end
    %Wczytanie wymuszen postepowych 
    no_wymp = str2num(fgetl(w));


    WymPostepowe = zeros(no_wymp,6);
    for i  = 1:no_wymp
        WymPostepowe(i,:) = str2num(fgetl(w));
        
    end
    %fprintf('Wymuszenia postępowe:\n')
    %disp(WymPostepowe)
    % Macierz masowa

    fclose(c);
    fclose(p);
    fclose(w);
    fclose(po);
end

%NewRaph.m
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
        disp('BŁĄD: Po założonej liczbie iteracji nie uzyskano zbieżności')
    end
end
%Predkosc.m
function [DQ] = Predkosc(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe,T_TEMP)
    %Funkcja rozwiązująca zadanie o prędkości
    
    FT = zeros(length(Q_TEMP), 1);
    Pozycja = 2*(size(ParyObrotowe, 1) + size(ParyPostepowe, 1))+1;
    
    %Wymuszenia postępowe
    for i=1:size(WymPostepowe,1)
        FT(Pozycja, 1) = -1*WymuszenieP(WymPostepowe(i, 2), 1, T_TEMP,ParyPostepowe, WymPostepowe);
        Pozycja = Pozycja + 1;
    end

    FQ = Jakobian(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe);
    DQ = -FQ \ FT;

end
%Przyspieszenie.m
function [D2Q] = Przyspieszenie(DQ_TEMP,Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe,T_TEMP)
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
        
        %Wzór 2.42 ze skryptu
        Gamma(Pozycja:Pozycja+1, 1) = Rot1 * sA * Dfi1^2 - Rot2 * sB * DFi2^2;
    
        Pozycja = Pozycja + 2;
    end
    
    %Pary postępowe
    for i=1:size(ParyPostepowe,1) 
        %Wzór 2.48     
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

        %Wzór 2.59
        Gamma(Pozycja, 1) = (Rot2 * V)'*(2 * Om * (DR2 - DR1) * DFi2 + (R2 - R1) * DFi2^2 - Rot1 * sA * (DFi2 - Dfi1)^2 ); 

        Pozycja = Pozycja + 1;
    end

    for i=1:size(WymPostepowe,1) %petla po wszystkich wymuszeniach postepowych    
        P1 = ParyPostepowe(WymPostepowe(i, 1), 1);
        P2 = ParyPostepowe(WymPostepowe(i, 1), 2);

        U = ParyPostepowe(WymPostepowe(i, 1), 4:5)';
        sA = ParyPostepowe(WymPostepowe(i, 1), 6:7)';

        [R1,~,Rot1] = Polozenia(Q_TEMP, P1);
        [R2,~,Rot2] = Polozenia(Q_TEMP, P2);

        [DR1,Dfi1,~]=Polozenia(DQ_TEMP,P1);
        [DR2,DFi2,~]=Polozenia(DQ_TEMP,P2);

        %Wzór 2.59
        Gamma(Pozycja, 1) = (Rot2 * U)'*(2 * Om * (DR2 - DR1) * DFi2 + (R2 - R1) * DFi2^2 - Rot1 * sA * (DFi2 - Dfi1)^2 ) +WymuszenieP(WymPostepowe(i, 2), 2, T_TEMP,ParyPostepowe, WymPostepowe); 

        Pozycja = Pozycja + 1;
    end

    %Obliczenie macierzy układu równań
    FQ = Jakobian(Q_TEMP, ParyObrotowe, ParyPostepowe, WymObrotowe, WymPostepowe);

    %Rozwiązanie
    D2Q = FQ \ Gamma;
end


%Wiezy.m
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
 
        %Wzór 2.18
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

        %Wzór 2.19
        F(Pozycja,1) = Fi1 - Fi2 - Fi0;
        
        %Wzór 2.22
        F(Pozycja+1,1) = (Rot2 * V)'*(R2 - R1 - Rot1 * sA) + V' * sB;
        
        Pozycja = Pozycja + 2;
    end

    %Stworzenie równań dla wymuszeń postępowych
    for i = 1:(size(WymPostepowe,1))
        [R1,~,Rot1] = Polozenia(Q_TEMP,ParyPostepowe(WymPostepowe(i,1),1));
        [R2,~,Rot2] = Polozenia(Q_TEMP,ParyPostepowe(WymPostepowe(i,1),2));
        U = ParyPostepowe(WymPostepowe(i,1),4:5)';
        sA = ParyPostepowe(WymPostepowe(i,1),6:7)';
        sB = ParyPostepowe(WymPostepowe(i,1),8:9)';
        
        % wzor 2.28 na wymuszenie w parze postepowej
        F(Pozycja,1) = (Rot2 * U)'*(R2 + Rot2 * sB - R1 - Rot1 * sA) - WymuszenieP(WymPostepowe(i,2),0,T_TEMP,ParyPostepowe, WymPostepowe);
        
        Pozycja = Pozycja + 1;
    end
end
%% jakobian

function [FQ] = Jakobian(Q_TEMP, ParyObrotowe, ParyPostepowe, ~, WymPostepowe)
    %Funkcja wyznaczająca Jakobian równań więzów 
    
    Om = [0 -1; 1 0];
    FQ = zeros(length(Q_TEMP),length(Q_TEMP));
    Pozycja = 1;

    %Pary obrotowe
    for i = 1:size(ParyObrotowe,1)
        sA = ParyObrotowe(i,3:4)';
        sB = ParyObrotowe(i,5:6)';
        [~,~,Rot1] = Polozenia(Q_TEMP,ParyObrotowe(i,1));
        [~,~,Rot2] = Polozenia(Q_TEMP,ParyObrotowe(i,2));
        
        %Je¿eli dany człon nie jest podstawą
        if ParyObrotowe(i,1) ~= 0 
            %Wzór 2.31
            FQ(Pozycja:Pozycja+1, ParyObrotowe(i,1)*3-2:ParyObrotowe(i,1)*3-1) = eye(2);
            
            %Wzór 2.32
            FQ(Pozycja:Pozycja+1,ParyObrotowe(i,1)*3) = Om * Rot1 * sA; 
        end
        
        %Jeżeli dany czlon nie jest podstawą
        if ParyObrotowe(i,2) ~= 0
            %Wzór 2.33
            FQ(Pozycja:Pozycja+1, ParyObrotowe(i,2)*3-2:ParyObrotowe(i,2)*3-1) = -eye(2); 
            
            %Wzór 2.34
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
        
        %Jeżeli dany człon nie jest podstawą
        if ParyPostepowe(i,1) ~= 0
            %Wzór 2.44
            FQ(Pozycja, ParyPostepowe(i,1)*3) = 1; 
            
            %Wzór 2.49
            FQ(Pozycja+1, ParyPostepowe(i,1)*3-2:ParyPostepowe(i,1)*3-1) = -(Rot2 * V)'; 
            
            %Wzór 2.50
            FQ(Pozycja+1, ParyPostepowe(i,1)*3) = -(Rot2 * V)' * Om * Rot1 * sA; 
        end
        
        %Je¿eli człon nie jest podstawą
        if ParyPostepowe(i,2) ~= 0
            %Wzór 2.46
            FQ(Pozycja, ParyPostepowe(i,2)*3) = -1; 
            
            %Wzór 2.51
            FQ(Pozycja+1, ParyPostepowe(i,2)*3-2:ParyPostepowe(i,2)*3-1) = (Rot2 * V)'; 
            
            %Wzór 2.52
            FQ(Pozycja+1, ParyPostepowe(i,2)*3) = -((Rot2 * V)' * Om * (R2 - R1 - Rot1 * sA)); 
        end
   
        Pozycja = Pozycja + 2;
    end

    %Wymuszenia postepowe
    for i = 1:size(WymPostepowe,1)
        U = ParyPostepowe(WymPostepowe(i,1),4:5)';
        sA = ParyPostepowe(WymPostepowe(i,1),6:7)';
        
        [R1,~,Rot1] = Polozenia(Q_TEMP,ParyPostepowe(WymPostepowe(i,1),1));
        [R2,~,Rot2] = Polozenia(Q_TEMP,ParyPostepowe(WymPostepowe(i,1),2));
        
        %Jeżeli człon nie jest podstawą
        if ParyPostepowe(WymPostepowe(i,1),1) ~= 0
            %Wzór 2.49
            FQ(Pozycja, ParyPostepowe(WymPostepowe(i,1),1)*3-2:ParyPostepowe(WymPostepowe(i,1),1)*3-1) = -(Rot2 * U)'; 
            
            %Wzór 2.50
            FQ(Pozycja, ParyPostepowe(WymPostepowe(i,1),1)*3) = -(Rot2 * U)' * Om * Rot1 * sA; % 
        end
        
        %Jeżeli człon nie jest podstawą
        if ParyPostepowe(WymPostepowe(i,1),2) ~= 0
            %Wzór 2.51
            FQ(Pozycja, ParyPostepowe(WymPostepowe(i,1),2)*3-2:ParyPostepowe(WymPostepowe(i,1),2)*3-1) = (Rot2 * U)'; 
            
            %Wzór 2.52
            FQ(Pozycja, ParyPostepowe(WymPostepowe(i,1),2)*3) = -(Rot2 * U)' * Om * (R2 - R1 - Rot1 * sA); 
        end
        Pozycja = Pozycja + 1;
    end
    
    %Sprawdzenie osobliwości Jakobianu
    if(det(FQ)==0)
        disp('Uwaga - Jakobian jest osobliwy')
    end
end

function [R, Fi, Roti] = Polozenia(Q,indeks)
    %Funkcja pobierająca aktualne współrzędne 

    if(indeks == 0) %Dla podstawy
        R = [0;0]; Fi = 0; Roti = Rot(Fi);
    else
        R = Q(3*indeks-2:3*indeks-1);
        Fi = Q(3*indeks);
        Roti = Rot(Fi);
    end
end
function delta = delta(Czlony,Punkty)
%%
    Punkty=Punkty(:);
    temp=Punkty(1:4:end);
    delta=zeros(length(Punkty)*3/4,1);
    for i=1:length(temp)
        delta(3*(i-1)+1:3*i)=[Punkty(4*i-2)-Czlony(3*temp(i)-2);Punkty(4*i-1)-Czlony(3*temp(i)-1);Punkty(4*i)-Czlony(3*temp(i))];
    end
    %%
end
function [pq,pdq,pd2q] = PunktPolozenie(Punkty,Czlony,Q,dQ,d2Q,nr_czlonu,nr_punktu)
%%
omega = [0 -1; 1 0];
    Punkty=Punkty(:);
    temp=Punkty(1:4:end);
    delta=zeros(length(Punkty)*3/4,1);
    for i=1:length(temp)
        delta(3*(i-1)+1:3*i)=[Punkty(4*i-2)-Czlony(3*temp(i)-2);Punkty(4*i-1)-Czlony(3*temp(i)-1);Punkty(4*i)-Czlony(3*temp(i))];
    end
pq=Q(3*nr_czlonu-2:3*nr_czlonu-1,:);
S_A=[delta(3*nr_punktu-2);delta(3*nr_punktu-1)];
%położenie
    for i=1:length(pq)
    pq(:,i)=pq(:,i)+Rot(Q(3*nr_czlonu,i))*S_A;
    end
%prędkość
    pdq=dQ(3*nr_czlonu-2:3*nr_czlonu-1,:);
for i=1:length(pdq)
    pdq(:,i)=pdq(:,i)+(omega*Rot(Q(3*nr_czlonu,i))*S_A*dQ(3*nr_czlonu,i));
end
%przyspieszenie
pd2q=d2Q(3*nr_czlonu-2:3*nr_czlonu-1,:);
for i=1:length(pd2q)
    pd2q(:,i)=pd2q(:,i)+(omega*Rot(Q(3*nr_czlonu,i))*S_A*d2Q(3*nr_czlonu,i))-(Rot(Q(3*nr_czlonu,i))*S_A*dQ(3*nr_czlonu,i)^2);
end
end
%%