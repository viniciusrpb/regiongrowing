function [kappa,smoothKappa,normKappa] = curvature(x,y)

    T = 25;
    dt = 0.1;
    tam = length(x);
    smoothKappa(1:tam+2,1:T) = 0;
    
    %Calcula angulo entre segmentos
    kappa(1:tam+2) = 0;
    xp(2:tam+1) = x;
    yp(2:tam+1) = y;
    xp(1) = x(end);
    yp(1) = y(end);
    xp(tam+2) = x(1);
    yp(tam+2) = y(1);

    %             figure,imshow(img,[]);
    %             hold on;
    %             plot(y(1),x(1),'--rs');
    %             pause;

    %Determina curvaturas
    for j = 2:tam+1
        alpha = (xp(j+1) - xp(j-1))/2;
        beta = yp(j+1) - 2*yp(j) + yp(j-1);
        epsn = (yp(j+1) - yp(j-1))/2;
        derta = xp(j+1) - 2*xp(j) + xp(j-1);
        kappa(j-1) = ((alpha*beta) - (epsn*derta)) / (( (alpha*alpha) + (epsn*epsn) + 1e-10 )^(3/2));
    end

    %Suaviza curva
    for j = 2:tam+1
        smoothKappa(j,1) = kappa(j-1);
        for t = 1:T-1
            bmais = smoothKappa(j+1,t)-smoothKappa(j,t);
            bmenos = smoothKappa(j,t) - smoothKappa(j-1,t);
            smoothKappa(j,t+1) = smoothKappa(j,t) + dt*( (g(bmais)*bmais) - (g(bmenos)*bmenos));
        end
    end

    NovoGamma(1:tam+2) = smoothKappa(1:tam+2,t);

    ff = double(1/tam);

    somaCurvs = sum(abs(NovoGamma));

    %Normaliza as curvaturas
    normKappa = NovoGamma ./ (ff * somaCurvs);
    
    kappa = kappa(2:end-1)';
    
    smoothKappa = smoothKappa(2:end-1,1:T);
    
    normKappa = normKappa(2:end-1)';
            
end

function y = g(x)
    y = 1/(1+(0.1*(x^2)));
end