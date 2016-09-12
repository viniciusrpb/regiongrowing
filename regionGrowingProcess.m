% Region growing algorithm based on Gonzalez, 2008: Digital Image Processing
% 3. ed.
% IMPORTANT: This implementation of region growing follows an iterative
% version proposed by Gonzalez et. al, 2008.

function bin = regionGrowingProcess(img,x,y,mu1,sigma1,mu2,sigma2,binTemp)

    % Take the weight and width
    [height,width] = size(binTemp);
    
    %Converts image to double
    f = double(img);
    
    %Compute the probabilities of each region and the seed
    fixed = 1/(2*sigma1*(pi^0.5));
    downfixed = 2.0*(sigma1*sigma1)+1e-40;
    sx = f(x,y);
    ps = fixed.*exp(-((sx-mu1).*(sx-mu1))./downfixed);
    pa = fixed.*exp(-((f-mu1).*(f-mu1))./downfixed);
    fixed = 1/(2*sigma2*(pi^0.5));
    downfixed = 2.0*(sigma2*sigma2)+1e-40;
    pb = fixed.*exp(-((f-mu2).*(f-mu2))./downfixed);
    
    % Detect the seed
    if numel(sx) == 1

        % Avoids division by zero
        r = pb./(pa+0.0000001);
        s = ps;
        
        % Normalize the threshold image for comparing probabilities
        maior = max(max(r));
        menor = min(min(r));
        rnorm = (r - menor)./(maior-menor);
        
        % Samples probabilities for each region using the pre-segmentation
        n1 = 1;
        n2 = 1;
        for alt = 1:height
            for larg = 1:width
                
                if(binTemp(alt,larg) ~= 0)
                    probRegion1(n1) = rnorm(alt,larg);
                    n1=n1+1;
                else
                    probRegion2(n2) =rnorm(alt,larg);
                    n2=n2+1;
                end
            end
        end
        
        % Computes the mean of probabilities values for each region
        M1 = mean(probRegion1);
        M2 = mean(probRegion2);
        
        % Obtain the satisfying regions that are able to receive similar
        % seeds
        avThreshold = (M1+M2)/2;
        
        % Keep the satisfying regions in image for the seed
        if M1 > M2
            si = (rnorm > avThreshold);
        else
            si = (rnorm < avThreshold);
        end

        s1 = s;
    end

    % Fill up ti with zero values
    ti=false(size(f));
    
    % For each seed
    for k = 1:length(s1)
        
        % Considers the current pixel denoted by k
        sv = s1(k);
        
        % Compute the dissimilarity between regions considering the seed
        s = power((f(:,:) - sv),2);
        ss = sqrt(double(s));
        
        % Normalize the dissimilarity image to [0,1]
        minss = min(min(ss));
        maxss = max(max(ss));
        normss = (ss - minss)./(maxss - minss);
        
        distancesInsidePreSeg = binTemp.*normss;
        dists = distancesInsidePreSeg(:);
        n=1;
        for l = 1:length(dists)
            if dists(l) > 0
                nonzeroDists(n) = dists(l);
                n=n+1;
            end
        end
        
        % Compute the mean value of
        %Include a tolerance to avoid zero values
        D0 = mean(nonzeroDists)+0.1;
        s = normss < min(D0,1.0);
        
        ti = ti | s;
    
    end
    
    %Reconstruct the image using using each seed label
    [g nr]=bwlabel(imreconstruct(si,ti,4));
    
    %Retrieve the regions for each seed
    ref = g(x,y);
    ti = (g == ref);
    
    bin = ti;
    
end
