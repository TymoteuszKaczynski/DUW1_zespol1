function plotp(T,dT,PunktPolozenie,type)
Time=0:dT:T;
if type=="x"
        figure()
        plot(Time,PunktPolozenie(1,:))
        xlabel('Czas [s]')
        ylabel('Położenie x punktu [m]')
        title("Położenie x punktu w czasie")
        grid on
elseif type=="y"
        figure()
        plot(Time,PunktPolozenie(2,:))
        ylabel('Położenie y punktu [m]')
        xlabel('Czas [s]')
        title('Położenie y punktu w czasie')
        grid on
end