function [descriptor,K] = cssDescriptor(x,y,sigma,numberOfPoints,thresh,dimensions)

    nroPontos = length(x);
    
    fator2 = 0.15;
    sigmaAtual(1) = 0.05;
    niveis = 1;
    while(sigmaAtual(niveis) < sigma)
        sigmaAtual(niveis+1) = double(sigmaAtual(niveis) + fator2);
        niveis = niveis + 1;
    end
    
    % Compute the CSS map
    [img,K] = cssMap(x,y,nroPontos,niveis,sigmaAtual,thresh);
    
    descriptor = constructDescriptor(uint8(img),numberOfPoints,sigmaAtual);

end

function descriptor = constructDescriptor(img,numberOfPoints,nivelSigma)

    [nroNiveis,nroPontos] = size(img);
    
    out = img(end:-1:1,:);
    
    maximas = 1;
    apagaFlag = 0;
    
    descriptor(1:numberOfPoints,1:2) = 0;
    
    for i = 4:nroNiveis-3
        for j = 4:nroPontos-3
            
            window = out(i,j-2:j+2);            
            
            if window(2) == 1 && window(3) == 1 && window(4) == 1
                descriptor(maximas,1) = j;
                descriptor(maximas,2) = nivelSigma(end-i+1);
                apagaFlag = 1;           
            else
                if window(2) == 1 && window(3) == 0 && window(4) == 1
                    descriptor(maximas,1) = j;
                    descriptor(maximas,2) = nivelSigma(end-i+1);
                    apagaFlag = 1;
                else
                    if window(2) == 1 && window(3) == 0 && window(4) == 0 && window(5) == 1
                        descriptor(maximas,1) = j;
                        descriptor(maximas,2) = nivelSigma(end-i+1);
                        apagaFlag = 1;
                    else
                        if window(1) == 1 && window(2) == 0 && window(3) == 0 && window(4) == 1
                            descriptor(maximas,1) = j;
                            descriptor(maximas,2) = nivelSigma(end-i+1);
                            apagaFlag = 1;
                        end
                    end
                end
            end

            if apagaFlag == 1
                maximas = maximas + 1;
                out(i-1:i+3,j-3:j+3) = 0;
                apagaFlag = 0;
            end
            
            if maximas == numberOfPoints+1
                break;
            end
              
        end
        
        if maximas == numberOfPoints+1
            break;
        end 
    end
    
    out = img(end:-1:1,:);
    
end

function [cssPlane,K] = cssMap(x,y,nroPontos,niveis,sigmaAtual,thresh)

    tamX = length(x);
    
    % Tamanho do filtro
    tam = tamX;
    if mod(tamX,2) == 0
        tam = tam + 1;
    end

    % Centralizacao do filtro
    inicio = floor(tam/2);

    g(1:tam) = 0;
    K(1:niveis,1:nroPontos) = 0;  
    
    cssPlane = zeros(niveis,nroPontos);
    
    % Processo de convolucao
    k = 1;
    
    
    
    for nivel = 1:niveis
        
        for c = -inicio:inicio
            g(c+inicio+1) = G(c,double(sigmaAtual(nivel)));
        end

        xs = conv1D(x,g);
        ys = conv1D(y,g);

        [kappa,kappaIterations,normKappa] = curvature(xs,ys);

        smoothKappa = kappa;


        % Encontra os cruzamentos em zero
        
        for j = 2:nroPontos-1
            
            if( smoothKappa(j-1)*smoothKappa(j+1) < 0 ) 
                
                abs1 = abs( smoothKappa(j-1) );
                abs2 = abs( smoothKappa(j+1) );
                
                %abs1 + abs2
                %pause;
                
                if (abs1 + abs2) > thresh
                    cssPlane(nivel,j) = double(1);
                    map(k,1) = j;
                    map(k,2) = sigmaAtual(nivel);
                    map(k,3) = smoothKappa(j);
                    k = k+1;
                end
            end
        end           
            
        
        K(nivel,:) = kappa;
        clear kappa;
        clear img;
        clear xs;
        clear ys;
        clear xxs;
        clear yys;
        
    end
   
end


function img = showSmoothImage(xs,ys,altura,largura)

    img = zeros(altura,largura);
    
    tam = length(xs);
    
    for i = 1:tam
        img(xs(i),ys(i)) = 255;
    end
    
end

function [xa,ya] = amostragemPontos(x,y)

    tam = length(x);

    space = ceil(tam/N);

    p = 1;
    for i = 1:space:tam
        xa(p) = x(i);
        ya(p) = y(i);
        p =p+1;
    end
    
    novoTam = length(xa);
    
    if novoTam < N    
        xa(end+1) = x(end);
        ya(end+1) = y(end);
    end

end

function GammaNorm = normaliza(Beta)

    tam= length(Beta);

    ff = double(1/tam);

    somaCurvs = sum(abs(Beta));

    %Normaliza as curvaturas
    GammaNorm = Beta ./ (ff * somaCurvs);

end

function y = G(x,sigma)

    esqB = sigma*sqrt(2*3.1416);
    esq = 1/esqB;
    x2 = - (x.*x);
    sig = 2*(sigma*sigma);
    intB = x2/sig;
    dir = exp(intB);
    y = esq*dir;

end

function y = conv1D(x,w)

    tamX = length(x);

    tamW = length(w);

    
    tiraCentro = tamW-1;
    
    inicio = tiraCentro/2 + 1;
    
    % Inicia vetor
    xx(1:tamX+tiraCentro) = 0;
    
    % Copia o vetor x para a parte adequada de xx
    xx(inicio:tamX+inicio-1) = x;
    
    % Replica o inicio
    xx(1:inicio-1) = x(end-inicio+2:end);
    
    % Replica o final
    xx(tamX+inicio:end) = x(1:inicio-1);
    
    % Copia x para a saida
    y = x;
    
    for i = inicio:tamX+inicio-1
        
        soma = 0;
        for j = -inicio+1:1:inicio-1
            soma = soma + w(inicio + j)*xx(i+j);
        end
        
        y(i-inicio+1) = soma;
    end
end