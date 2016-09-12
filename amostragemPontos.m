function [xa,ya] = amostragemPontos(x,y,N)

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
    
    xa = xa';
    ya = ya';

end