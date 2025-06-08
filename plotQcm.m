function plotQcm(T,dT,Q,num,type)
Time=0:dT:T;
if type=="x"
        figure()
        plot(Time,Q(3*num-2,:))
        xlabel('Czas [s]')
        ylabel('Położenie x punktu [m]')
        title("Położenie x punktu w czasie")
        grid on
elseif type=="y"
        figure()
        plot(Time,Q(3*num-1,:))
        ylabel('Położenie y punktu [m]')
        xlabel('Czas [s]')
        title('Położenie y punktu w czasie')
        grid on

elseif type=="fi"
        figure()
        plot(Time,Q(3*num,:))
        ylabel('Położenie fi punktu [rad]')
        xlabel('Czas [s]')
        title('Położenie fi punktu w czasie')
        grid on
end
end