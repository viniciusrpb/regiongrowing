% Implementation of the Curvature Scale Space Descriptor
% 
% INPUT: binary image img
% OUTPUT: the threshold cut at the CSS Map

function cutHere = cssSmoothing(img)

    % Create background image
    area = sum(sum(img));
    
    boundaries = Extracao_Contorno(img);
    x = boundaries(1,:);
    y = boundaries(2,:);
    
    x = x';
    y = y';
    
    qtdeAmostragem = 200;
    [xa,ya] = amostragemPontos(x,y,qtdeAmostragem);
    [cssMap,cssK]  = cssDescriptor(xa,ya,32.0,10,0.0,area);
    
    valid = cssMap(:,2) > 7;
    cutHere = sum(valid);


end