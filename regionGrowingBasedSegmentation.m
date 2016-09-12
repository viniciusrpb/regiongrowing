% INSTITUTO DE CIENCIAS MATEMATICAS E DE COMPUTACAO
% UNIVERSITY OF SAO PAULO, SAO CARLOS
% MATLAB CODE
% IMAGE SEGMEN
% LAST UPDATE: May, 12th 2016

function regionGrowingBasedSegmentation(path,pathSaida)

    pngPath = sprintf('%s*.png',path);
    list = dir(pngPath);
    listSize = length(list);
    
    %Flag write intermediate images to disk
    flag = 1;
    
    for i =1:listSize
        
        % Read image
        pathh = sprintf('%s%s',path,list(i).name);
        fprintf('(%d) - IMAGE: %s\n',i,pathh);
        
        % Transform to double
        imgRaw = double(imread(pathh));
        
        fprintf('Anisotropic diffusion: ');
        
        % Smooth and normalize the RGB image
        imgF = filtroEDP(uint8(imgRaw),0.01,0.2,20,0.01);
        
        fprintf('OK \n\n');
        
        
        [height,width,channels] = size(imgRaw);

        imgRGB(1:height,1:width,1:3) = 0;

        imgF = double(imgF);
        for c = 1:3
            maior = max(max(imgF(:,:,c)));
            menor = min(min(imgF(:,:,c)));
            imgRGB(:,:,c) = ( imgF(2:end-1,2:end-1,c) - menor) ./ (maior - menor + 1e-10);
        end

        % Eigenvalues image computation
        fprintf('Pre-segmentation image: ');
        imgRGB = double(imgRGB);
        [color_eigenValuesImageRaw,color_eigenVectorImage,color_values,color_vectors]= localPrincipalComponents(imgRGB,7,2,'simplex1testes/');
        for c = 1:3
            maior = max(max(color_eigenValuesImageRaw(:,:,c)));
            menor = min(min(color_eigenValuesImageRaw(:,:,c)));
            color_eigenValuesImage(:,:,c) = ( color_eigenValuesImageRaw(:,:,c) - menor) ./ (maior - menor + 1e-10);
        end
        

        % % % THRESHOLD OF THE EIGENVALUE IMAGE
        thresholdEigenValue = mean(mean(color_eigenValuesImage(:,:,3)));
        
        % Seeds computation and placement
        binMask = color_eigenValuesImage(:,:,3) > thresholdEigenValue;
        
        % Erosion procedure to the Mask Image to remove holes
        initialMask = imerode(binMask,strel('disk',3));
        
        % Prepare the mask
        binTemp = zeros(size(initialMask,1),size(initialMask,2));
        binTemp(10:size(initialMask,1)-10,10:size(initialMask,2)-10) = initialMask(10:size(initialMask,1)-10,10:size(initialMask,2)-10);
        fprintf('OK \n\n');
        
        fprintf('Determine the seeds: ');
        
        % Find the seeds for region growing
        [localInfo,perReference] = findSeeds(binTemp);    
        
        fprintf('OK \n\n');
        
        fprintf('Contrast enhancement of hue channel: ');
        % The enhancement process
        referenceHue =enhacementProcess(imgRGB,perReference);    
        fprintf('OK \n\n');
        
        % Show the number of found seeds
        numberOfSeeds = length(localInfo);
        
        % Images for illustration
        finalImg = zeros(height,width);
        seedsSubareas= zeros(height,width);
        outputSeeds(:,:,1:3) = uint8(imgRaw);

        % Samplind procedure ilustration
        imgSample = imgRaw;
        maxN = round(height*width*0.1);
        
        fprintf('Sampling procedure: ');
        
        nb=1;
        n=1;
        for a = 1:height
            for b=1:width
                if binTemp(a,b) == 1
                    imgSample(a,b,1) = 255;
                    tempVector(n) = double(referenceHue(a,b));
                     n=n+1;
                else
                    if nb < maxN
                        imgSample(a,b,2) = 255;
                        background(nb) = double(referenceHue(a,b));
                        nb=nb+1;
                    end
                end
                nb=nb+1;
            end
        end
        fprintf('OK \n\n');

        % For each seed
        for seed = 1:numberOfSeeds

            maskImg = localInfo{seed}.area;
            centroids(seed,1) = localInfo{seed}.cx;
            centroids(seed,2) = localInfo{seed}.cy;
            
            outputSeeds(centroids(seed,2)-6:centroids(seed,2)+6,centroids(seed,1)-6:centroids(seed,1)+6,1) = 255;
            outputSeeds(centroids(seed,2)-6:centroids(seed,2)+6,centroids(seed,1)-6:centroids(seed,1)+6,2) = 0;
            outputSeeds(centroids(seed,2)-6:centroids(seed,2)+6,centroids(seed,1)-6:centroids(seed,1)+6,3) = 0;

            %Compute the mean and standard deviation
            algaeMean = mean(tempVector);
            algaestandardDeviation = std(tempVector);
            backgroundMean = mean(background);
            backgroundstandardDeviation = std(background);
            
            % Plot Gaussians
            if flag == 1
                x = -300:.5:300;
                norm = normpdf(x,algaeMean,algaestandardDeviation);
                norm2 = normpdf(x,backgroundMean,backgroundstandardDeviation);
                figure,plot(x,norm,'r',x,norm2,'g','LineWidth',3);
            end

            fprintf('Region growing: ');
            % The region growing process
            regrow = regionGrowingProcess(referenceHue,round(centroids(seed,2)),round(centroids(seed,1)),algaeMean,algaestandardDeviation,backgroundMean,backgroundstandardDeviation,maskImg);

            fprintf('OK \n\n');
            
            % Get perimeter information
            allBoundaries = bwboundaries(regrow,'noholes');

            x = allBoundaries{1}(:,1);

            tamx = length(x);

            perimeter = tamx;
            
            % The rolling ball transformation
            iterationsNumber = max(1,cssSmoothing(regrow));
            if perimeter > 250 && iterationsNumber <= 3 
                radius = 6;
            else
                radius = 2;
            end
            dilImg = imclose(regrow,strel('disk',radius));
            finalSegmRolBal = double(dilImg > 0);
            
            allBoundaries = bwboundaries(finalSegmRolBal,'noholes');

            x = allBoundaries{1}(:,1);
            y = allBoundaries{1}(:,2);

            tamx = length(x);

            for l = 1:tamx
                finalImg(x(l),y(l)) = 1;
            end

            finalImg = imfill(finalImg,'holes');

            close all;
            clear finalImg2;
            clear finalImg3;
            clear regionOfInterest;
            clear regionOfInterest2;
            clear binFinal;
            clear bin2;
            clear imgN;
            clear valueDist2;
            clear roi_coord;
            clear maskImg;
        end
        
        if flag ==1
            pathf = sprintf('%s%s_seeds.png',pathSaida,list(i).name(1:end-4));
            imwrite(outputSeeds,pathf);

            pathf = sprintf('%s%s_segmentation.png',pathSaida,list(i).name(1:end-4));
            imwrite(uint8(finalImg*255),pathf);
        end
    end
       
end

% Compute the seeds making sure they are inside algae regions
function [localInfo,maiorDos] = findSeeds(referenceImg)

    initialImage = referenceImg;
    [height,width] = size(initialImage);

    regionBoundaries = bwboundaries(initialImage,'noholes');

    for j = 1:length(regionBoundaries)

        x = regionBoundaries{j}(:,1);
        y = regionBoundaries{j}(:,2);

        tamx = length(x);

        boundaries = zeros(height,width);
        for l = 1:tamx

            boundaries(x(l),y(l)) = 1;

        end

        area = imfill(boundaries,'holes');
        stats = regionprops(area, 'Perimeter');
        perimeters(j) = stats.Perimeter;

    end

    candidates = perimeters > 150
    maiorDos = max(perimeters);
    [sortedPerimeters,indexes] = sort(perimeters.*candidates,'descend');

    realRegions = sum(sum(sortedPerimeters > 0));

    % Set the disk radius for the erosion operation
    numberOfRegions = sum(realRegions);
    seedsNumber = numberOfRegions;

    % Compute the centroid for each region
    ss=1;
    for s = 1:seedsNumber

        x = regionBoundaries{indexes(s)}(:,1);
        y = regionBoundaries{indexes(s)}(:,2);

        tamx = length(x);

        fronteiras = zeros(height,width);
        for l = 1:tamx
            fronteiras(x(l),y(l)) = 1;
        end

        area = imfill(fronteiras, 'holes');

        elem = strel('disk',1);
        er = area;
        soma = sum(sum(er));
        thick = 0;
        while(soma > 0)
            er = imerode(er,elem);
            soma = sum(sum(er));
            thick = thick + 1;
        end

        avgDisk = ceil(thick/2);
        
        subarea = imerode(area,strel('disk',avgDisk));

        if(sum(sum(subarea))) > 0

            todasBordas2 = bwboundaries(subarea);

            x = todasBordas2{1}(:,1);
            y = todasBordas2{1}(:,2);

            [value,ind] = min(x);
            cx = y(ind);
            cy = x(ind);


            localInfo{ss}.cx = cx;
            localInfo{ss}.cy = cy;
            localInfo{ss}.area = area;
            ss=ss+1;
        end
        
        clear x;
        clear y;
        clear todasBordas2;
        clear fronteiras;
        clear subarea;
        clear roi;
        clear area;

    end
end
