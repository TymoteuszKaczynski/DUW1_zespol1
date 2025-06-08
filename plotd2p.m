function plotd2p(T,dT,PunktPrzyspieszenie,type)
Time=0:dT:T;
if type=="x"
        figure()
        plot(Time,PunktPrzyspieszenie(1,:))
        xlabel('Czas [s]')
        ylabel('Przyspieszenie x punktu [m/s^2]')
        title("Przyspieszenie x punktu w czasie")
        grid on
elseif type=="y"
        figure()
        plot(Time,PunktPrzyspieszenie(2,:))
        ylabel('Przyspieszenie y punktu [m/s^2]')
        xlabel('Czas [s]')
        title('Przyspieszenie y punktu w czasie')
        grid on
end