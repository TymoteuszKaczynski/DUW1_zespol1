function plotDQcm(T,dT,DQ,num,type)
Time=0:dT:T;
if type=="x"
        figure()
        plot(Time,DQ(3*num-2,:))
        xlabel('Czas [s]')
        ylabel('Prędkość x punktu [m/s]')
        title("Prędkość x punktu w czasie")
        grid on
elseif type=="y"
        figure()
        plot(Time,DQ(3*num-1,:))
        ylabel('Prędkość y punktu [m/s]')
        xlabel('Czas [s]')
        title('Prędkość y punktu w czasie')
        grid on

elseif type=="fi"
        figure()
        plot(Time,DQ(3*num,:))
        ylabel('Prędkość fi punktu [rad/s]')
        xlabel('Czas [s]')
        title('Prędkość fi punktu w czasie')
        grid on
end
end