function plotdp(T,dT,PunktPredkosc,type)
Time=0:dT:T;
if type=="x"
        figure()
        plot(Time,PunktPredkosc(1,:))
        xlabel('Czas [s]')
        ylabel('Prędkość x punktu [m/s]')
        title("Prędkość x punktu w czasie")
        grid on
elseif type=="y"
        figure()
        plot(Time,PunktPredkosc(2,:))
        ylabel('Prędkość y punktu [m/s]')
        xlabel('Czas [s]')
        title('Prędkość y punktu w czasie')
        grid on
end