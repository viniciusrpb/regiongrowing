%Filtro de Suavizacao e Eliminacao de Ruidos baseado em Equacoes
% Diferenciais Parciais, proposto pela prof. Dra. Celia Aparecida
% Zorzo Barcelos.
%Implementado por Vinicius Ruela Pereira Borges, 08/01/2009

function U = filtroEDP(img,k,dt,T,LAMBDA)
    % k = fator de suavizacao (quanto menor, maior eh a suavizacao)
    % dt = ponderacao do filtro [0.1 0.4] - estabilidade
    % T = quantidade de iteracoes no tempo
    %LAMBDA = pondera o balanceamento (segundo termo da EDP)
    %         valor esperado entre [0 1]

    [altura,largura,channels] = size(img);
    
    if channels == 1

        %Inicia matriz de saida
        Un_1 = zeros(altura+4,largura+4);

        %Ajusta tamanho da imagem
        I = padarray(img,[2 2],'replicate');

        %Condicoes Inicial e Contorno
        U = I;

        %T iteracoes no tempo
        for t = 1:T

            %Percorre a imagem inteira
            for j = 3:altura+2
                for i = 3:largura+2

                    %Aproximacao dos Termos por Diferencas Finitas
                    ux = (U(j,i+1) - U(j,i-1))/2;
                    uy = (U(j+1,i) - U(j-1,i))/2;
                    uxx = (U(j,i+1) - 2*U(j,i) + U(j,i-1));
                    uyy = (U(j+1,i) - 2*U(j,i) + U(j-1,i));
                    uxy = (U(j+1,i+1) - U(j-1,i+1) - U(j+1,i-1) + U(j-1,i-1))/4;

                    %Aproximacao dos Termos por Derivadas de Segunda Ordem por
                    %Linha
                    dudx = abs((U(j-1,i+1) + (2*U(j,i+1)) + U(j+1,i+1)) - (U(j-1,i-1) + (2*U(j,i-1)) + U(j+1,i-1)));
                    dudy = abs((U(j+1,i-1) + (2*U(j+1,i)) + U(j+1,i+1)) - (U(j-1,i-1) + (2*U(j-1,i)) + U(j-1,i+1)));


                    %Calculo de g
                    s = ux^2 + uy^2 + (1e-10);
                    g = 1/(1+(k*(dudx +dudy)));

                    u2x = ux*ux;
                    u2y = uy*uy;

                    %Calculo do TAL, onde 0.000001 evita divisao por zero
                    TAL = g*(((u2x*uyy) - (2*(ux*(uy*uxy))) + (u2y*uxx))/(u2x + u2y + (1e-10))) + LAMBDA*(1-g)*(U(j,i) - I(j,i));

                    %Calculo da variavel temporal no tempo n+1
                    Un_1(j,i) = U(j,i) + dt*TAL;

                end
            end
            %Passa o tempo
            U = Un_1;
        end

        %Retoma �rea inicial da imagem
        U = Un_1(2:altura+3,2:largura+3);
    
    else
        
        %Inicia matriz de saida
        for c = 1:channels
            Un_1(:,:,c) = zeros(altura+4,largura+4);
        end
        
        %Ajusta tamanho da imagem
        I = padarray(double(img),[2,2],'replicate');

        
        %Para cada canal da imagem
        for c = 1:channels
            
            %Condicoes Inicial e Contorno
            U = I(:,:,c);
        
            %T iteracoes no tempo
            for t = 1:T
                
                %Percorre a imagem inteira
                for j = 3:altura+2
                    for i = 3:largura+2

                        %Aproximacao dos Termos por Diferencas Finitas
                        ux = (U(j,i+1) - U(j,i-1))/2;
                        uy = (U(j+1,i) - U(j-1,i))/2;
                        uxx = (U(j,i+1) - 2*U(j,i) + U(j,i-1));
                        uyy = (U(j+1,i) - 2*U(j,i) + U(j-1,i));
                        uxy = (U(j+1,i+1) - U(j-1,i+1) - U(j+1,i-1) + U(j-1,i-1))/4;

                        %Aproximacao dos Termos por Derivadas de Segunda Ordem por
                        %Linha
                        dudx = abs((U(j-1,i+1) + (2*U(j,i+1)) + U(j+1,i+1)) - (U(j-1,i-1) + (2*U(j,i-1)) + U(j+1,i-1)));
                        dudy = abs((U(j+1,i-1) + (2*U(j+1,i)) + U(j+1,i+1)) - (U(j-1,i-1) + (2*U(j-1,i)) + U(j-1,i+1)));


                        %Calculo de g
                        %s = ux^2 + uy^2 + (1e-10);
                        g = 1/(1+(k*(dudx +dudy)));

                        u2x = ux*ux;
                        u2y = uy*uy;

                        %Calculo do TAL, onde 0.000001 evita divisao por zero
                        TAL = g*(((u2x*uyy) - (2*(ux*(uy*uxy))) + (u2y*uxx))/(u2x + u2y + (1e-10))) + LAMBDA*(1-g)*(U(j,i) - I(j,i,c));

                        %Calculo da variavel temporal no tempo n+1
                        Un_1(j,i) = U(j,i) + dt*TAL;

                    end
                end
                %Passa o tempo
                U = Un_1;
            end
            U_final(:,:,c) = U(:,:);
        end

        %Retoma area inicial da imagem
        clear U;
        U = U_final(2:altura+3,2:largura+3,1:channels);
        
        
    end
%     %Ilustra resultados
%     figure,imshow(img,[]);
%     title('Imagem Original');
% 
%     %Imagem Filtrada com EDP
%     figure,imshow(U,[]);
%     title('Imagem Filtrada');
    
end

