function [eigenValuesImage,eigenVectorImage,eigenvaluesorig,eigenvectorsorig]= localPrincipalComponents(img,radius,casee,path)

    [height,width,channels] = size(img);
    
    middle = floor(radius/2);
    
    eigenValuesImage = zeros(height,width);
    eigenVectorImage = zeros(height,width);
    
    if channels == 1
        
        if casee == 1
    
            places(1:2,1:2) = 0;

            %V = radius*radius;

            no_dims=2;

            imgpad = padarray(double(img),[middle,middle],'replicate');

            for x=middle+1:height-middle
                for y=middle+1:width-middle

                    n = 1;
                    for rx = -middle:1:middle
                        for ry = -middle:1:middle

                            xx = x+rx;
                            yy = y+ry;

                            if imgpad(xx,yy) ~= 0
                                places(n,1) = xx/middle;
                                places(n,2) = yy/middle;
                                n=n+1;
                            end

                        end
                    end
                    n=n-1;

                    M(1) = mean(places(:,1));
                    M(2) = mean(places(:,2));

                    soma(1:2,1:2) = 0;

                    for i = 1:n
                        soma = soma + ( (  (double(places(i,:) - M )  )' * ( double(places(i,:)) - M ) ));
                    end

                    Cov(:,:) = (1/(n-1))*soma;

                    Cov(isnan(Cov)) = 0;
                    Cov(isinf(Cov)) = 0;
                    [eigenvectorsorig, eigenvaluesorig] = eig(Cov);
                    [eigenvalues,ind] = sort(diag(eigenvaluesorig),'descend');
                    eigenvectors = eigenvectorsorig(:,ind(1:no_dims));
                    eigenvalues = eigenvalues(1:no_dims);

                    %pc = eigenvectors(:,ind(1:no_dims)); 

                    eigenValuesImage(x-middle,y-middle,1) = eigenvalues(1);
                    eigenVectorImage(x-middle,y-middle,1) = eigenvectors(1);
                    eigenValuesImage(x-middle,y-middle,2) = eigenvalues(2);
                    eigenVectorImage(x-middle,y-middle,2) = eigenvectors(2);
                    %localAverage(x-middle,y-middle,1:channels) = sumLocal(1:channels)./N;

                end
            end
        else
            
            places(1:2) = 0;

            %V = radius*radius;

            no_dims=1;

            imgpad = padarray(double(img),[middle,middle],'replicate');

            for x=middle+1:height-middle
                for y=middle+1:width-middle

                    n = 1;
                    for rx = -middle:1:middle
                        for ry = -middle:1:middle

                            xx = x+rx;
                            yy = y+ry;

                            %if imgpad(xx,yy) ~= 0
                                places(n) = imgpad(xx,yy);
                                %places(n,2) = yy/middle;
                                n=n+1;
                            %end

                        end
                    end
                    n=n-1;

                    M = mean(places);
                    %M(2) = mean(places(:,2));

                    soma(1:2,1:2) = 0;

                    for i = 1:n
                        soma = soma + ( (  (double(places(i) - M )  )' * ( double(places(i)) - M ) ));
                    end

                    Cov(:,:) = (1/(n-1))*soma;

                    Cov(isnan(Cov)) = 0;
                    Cov(isinf(Cov)) = 0;
                    [eigenvectorsorig, eigenvaluesorig] = eig(Cov);
                    [eigenvalues,ind] = sort(diag(eigenvaluesorig),'descend');
                    eigenvectors = eigenvectorsorig(:,ind(1:no_dims));
                    eigenvalues = eigenvalues(1:no_dims);

                    %pc = eigenvectors(:,ind(1:no_dims)); 

                    eigenValuesImage(x-middle,y-middle,1) = eigenvalues(1);
                    eigenVectorImage(x-middle,y-middle,1) = eigenvectors(1);
                    %localAverage(x-middle,y-middle,1:channels) = sumLocal(1:channels)./N;

                end
            end
        end
    else
        
        places(1:2,1:3) = 0;

        %V = radius*radius;

        no_dims=3;

        imgpad = padarray(double(img),[middle,middle],'replicate');

        for x=middle+1:height-middle
            for y=middle+1:width-middle

                n = 1;
                for rx = -middle:1:middle
                    for ry = -middle:1:middle

                        xx = x+rx;
                        yy = y+ry;

                        %if imgpad(xx,yy) ~= 0
                            places(n,1) = imgpad(xx,yy,1);%xx/middle;
                            places(n,2) = imgpad(xx,yy,2);
                            places(n,3) = imgpad(xx,yy,3);
                            n=n+1;
                        %end

                    end
                end
                n=n-1;

                M(1) = mean(places(:,1));
                M(2) = mean(places(:,2));
                M(3) = mean(places(:,3));
                
      

                soma(1:3,1:3) = 0;
                
 

                for i = 1:n
                    soma = soma + ( (  (double(places(i,:) - M )  )' * ( double(places(i,:)) - M ) ));
                end

                Cov(:,:) = (1/(n-1))*soma;

                Cov(isnan(Cov)) = 0;
                Cov(isinf(Cov)) = 0;
                [eigenvectorsorig, eigenvaluesorig] = eig(Cov);
                [eigenvalues,ind] = sort(diag(eigenvaluesorig),'descend');
                eigenvectors = eigenvectorsorig(:,ind(1:no_dims));
                eigenvalues = eigenvalues(1:no_dims);

                %pc = eigenvectors(:,ind(1:no_dims)); 

                eigenValuesImage(x-middle,y-middle,1) = eigenvalues(1);
                eigenVectorImage(x-middle,y-middle,1) = eigenvectors(1);
                eigenValuesImage(x-middle,y-middle,2) = eigenvalues(2);
                eigenVectorImage(x-middle,y-middle,2) = eigenvectors(2);
                eigenValuesImage(x-middle,y-middle,3) = eigenvalues(3);
                eigenVectorImage(x-middle,y-middle,3) = eigenvectors(3);
                %localAverage(x-middle,y-middle,1:channels) = sumLocal(1:channels)./N;

            end
        end
        
    end
    
%     figure,imshow(eigenVectorImage,[]);
%     title('Eigen vector');
%     figure,imshow(eigenValuesImage,[]);
%     title('Eigen values');
%     pause;
%     ma = max(max(eigenValuesImage));
%     mi = min(min(eigenValuesImage));
%     eigVa = (eigenValuesImage - mi)./(ma-mi);
%     paths = sprintf('eigenimages/%s_eigvalue.png',path);
%     imwrite(uint8(eigVa*255),paths);
%     
%     ma = max(max(eigenVectorImage));
%     mi = min(min(eigenVectorImage));
%     eigVe = (eigenVectorImage - mi)./(ma-mi);
%     paths = sprintf('eigenimages/%s_eigvector.png',path);
%     imwrite(uint8(eigVe*255),paths);
    close all;
    
end