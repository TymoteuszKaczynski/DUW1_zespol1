function plotD2Qcm(T,dT,D2Q,num,type)
Time=0:dT:T;
if type=="x"
        figure()
        plot(Time,D2Q(3*num-2,:))
        xlabel('Czas [s]')
        ylabel('Przyspieszenie x punktu [m/s^2]')
        title("Przyspieszenie x punktu w czasie")
        grid on
elseif type=="y"
        figure()
        plot(Time,D2Q(3*num-1,:))
        ylabel('Przyspieszenie y punktu [m/s^2]')
        xlabel('Czas [s]')
        title('Przyspieszenie y punktu w czasie')
        grid on

elseif type=="fi"
        figure()
        plot(Time,D2Q(3*num,:))
        ylabel('Przypsieszenie fi punktu [rad/s^2]')
        xlabel('Czas [s]')
        title('Przyspieszenie fi punktu w czasie')
        grid on
end
end