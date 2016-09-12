%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Algoritmo de extração do contorno de uma imagem binarizada          %%
%%         **********************************************************          %%
%%       O método aqui implementado extrai o contorno do objeto na imagem,     %%
%%       conforme algoritmo de rastreamento de fronteira de Moore proposto     %%
%%       no livro do Gonzalez (pág. 524).                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Entrada:                                                                    %%
%%   - IMGP:  imagem binarizada com um objeto.                                 %%
%% Saída:                                                                      %%
%%   - CONTORNO: coordenadas de contorno da imagem.                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CONTORNO = Extracao_Contorno(IMGP)

        DIM = size(IMGP); %% Obtem dimensão da imagem.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Geração do contorno pelo método proposto em Gonzalez.                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Criando imagem com Padding
        IMGP = [zeros(1,DIM(2)); IMGP; zeros(1,DIM(2))];
        IMGP = [zeros(DIM(1)+2,1), IMGP, zeros(DIM(1)+2,1)];
        DIMP = size(IMGP);

        % Máscara definida a partir da variaçao do pixel de contorno atual
        MASC = [ 0,-1,-1,-1, 0, 1, 1, 1 ;
                -1,-1, 0, 1, 1, 1, 0,-1 ];

        %% Busca o 1o pixel do objeto (pixel + acima e + a esquerda)
        for ix=1:DIMP(1)      % ix e iy sao as coordenadas da imagem
            for iy=1:DIMP(2)
                if IMGP(ix,iy) == 1
                   break;
                end
            end
            if iy < DIMP(2)
               break;
            end
        end
        %% Configura os parametros iniciais do objeto
        CONTORNO = [ix;iy]; %% Guarda o 1o pixel
        B = CONTORNO;
        C = 0; % Posicao da mascara imediatamente anterior ao pixel

        %% Aplicando o algoritmo do ceginho
        contador = 0; % recurso para evitar o loop eterno
        while ((ix < DIMP(1)) || (iy < DIMP(2))) && (contador <= size(CONTORNO,2))
              %% Pecorre a mascara a partir de C até encontrar um novo pixel
              for IND=C:(C+7)
                  POS = mod(IND,8)+1; % Recurso para garantir o pecorrimento ciclico
                  if IMGP(B(1)+MASC(1,POS), B(2)+MASC(2,POS)) == 1
                      break;
                  end
              end
              if IND == (C+8) % O ceguinho deu uma volta completa sem achar ninguém
                  %% Busca novo 1o pixel do objeto
                  for ix=B(1):DIMP(1)
                      if ix == B(1)
                          INI = B(2)+1;
                      else
                          INI = 1;
                      end
                      for iy=INI:DIMP(2)
                          if IMGP(ix,iy) == 1
                              break;
                          end
                      end
                      if iy < DIMP(2)
                          break;
                      end
                  end
                  %% Configura os parametros iniciais do objeto
                  CONTORNO = [ix;iy]; %% Guarda o 1o pixel
                  B = CONTORNO;
                  C = 0; % Posicao da mascara imediatamente anterior ao pixel
                  continue;
              end
              %% Busca por repetições no arquivo de contorno
              achei = 0;
              for fx = 1:size(CONTORNO,2)
                  if ((CONTORNO(1,fx) == B(1)+MASC(1,POS)) && (CONTORNO(2,fx) == B(2)+MASC(2,POS)))
                      achei = 1;
                      break;
                  end
              end
              if achei == 0 %% Se nao achar
                  CONTORNO = [CONTORNO , [B(1)+MASC(1,POS); B(2)+MASC(2,POS)]];
                  contador = 0;
              else
                  if (fx == 2) && ((CONTORNO(1,1) == B(1)) && (CONTORNO(2,1) == B(2))) % Criterio de parada (2 ultimos pontos = 2 primeiros)
                      break;
                  else
                      contador = contador+1;
                  end
              end
              B = [B(1)+MASC(1,POS); B(2)+MASC(2,POS)];
              C = mod((POS+5),8)-(1-mod(POS,2));
        end
        %% Subtrai a linha e coluna incluída pelo padding
        CONTORNO = CONTORNO-1;
        %% ordena as coordenadas do contorno
        %CONTORNO = sortrows(CONTORNO');
        %CONTORNO = CONTORNO';            

        %% Limpa variaveis auxiliares
        clear C B fx achei IND MASC POS IMGP DIM DIMP ix iy contador;
