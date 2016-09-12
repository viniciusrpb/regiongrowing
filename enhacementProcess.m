function Ien = enhacementProcess(smoothedRGB,perimeter)
    
    % Gets the smoothed HSV model
    hsvSmoothed = rgb2hsv(double(smoothedRGB));
    
    % Gets the values channel
    valueChannel =hsvSmoothed(:,:,3);
    
    % Gets the hue channel
    hue = hsvSmoothed(:,:,1);
    
    % Computes the equalized image to 64 gray intensities
    equalizedImg = histeq(uint8(valueChannel*255),64);

    signal = 1;
    if perimeter < 250
        signal = 0;
    end
    
    if perimeter < 500
        threshold = 1;
    else
        if perimeter < 1000
            threshold = 5;
        else
            if perimeter < 1500
                threshold = 9;
            else
                threshold = 12;
            end
        end
    end
    
    % Normalization to [0,1]
    maior = max(max(hue));
    menor = min(min(hue));
    huen = ( hue - menor) ./ (maior - menor+1e-20);

    % Obtain the regions to be enhanced
    BEQ = equalizedImg >= threshold;
    
    % Enhance the hue channel
    if signal > 0
        hh = huen + 2*huen.*BEQ;
    else
        hh = huen+BEQ;
    end

    % Normalization to [0,1]
    higher = max(max(hh));
    fewer = min(min(hh));
    enhancedNorm = ( hh - fewer) ./ ( higher - fewer);
    
    Ien = uint8(enhancedNorm*255);

end