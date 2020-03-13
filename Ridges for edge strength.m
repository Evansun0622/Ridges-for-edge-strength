[IA, std_dev, IB] = problemOne(imread('Jayce_2.jpg'))
[xx, xy, yy, eigenvalues, image_e] = problemThree(IB)
function [xx, xy, yy, eigenvalues, image_e] = problemThree(blurred_diagonal)
    [Amplitude, Phase] = AmpPhaseDFT(blurred_diagonal);
    for j = 1:65
        for k = 1:128
            if k < 66
                F(j,k) = Amplitude(j,k)*exp(-((k-1)/128)^2*(2*pi)^2);
            else
                F(j,k) = Amplitude(j,k)*exp(-((k-129)/128)^2*(2*pi)^2);
            end
        end
    end

    for k = 1:128
        for j = 1:65
            F(j,k) = F(j,k)*exp(-((j-1)/128)^2*(2*pi)^2);
        end
    end
    
    for j = 1:65
        for k = 1:128
            if k < 66
                xx(j,k) = F(j,k)*(-((k-1)/128)^2);
            else
                xx(j,k) = F(j,k)*(-((k-129)/128)^2);
            end
        end
    end
    xx(1,1) = 0;
    xx(65,1) = 0;
    xx(1,65) = -0.25 * F(1,65);
    xx(65,65) = -0.25 * F(65,65);
    xx = ReconfromAmpPhase(xx,Phase);
    for j = 1:65
        for k = 1:128
            if k < 66
                xy(j,k) = F(j,k)*(-(k-1)/128)*((j-1)/128);
            else
                xy(j,k) = F(j,k)*(-(k-129)/128)*((j-1)/128);
            end
        end
    end
    xy(1,1) = 0;
    xy(65,1) = 0;
    xy(1,65) = 0;
    xy(65,65) = -0.25 * F(65,65);
    xy = ReconfromAmpPhase(xy,Phase);
    
    for j = 1:65
        for k = 1:128
            yy(j,k) = F(j,k)*(-((j-1)/128)^2);
        end
    end
    yy(1,1) = 0;
    yy(1,65) = 0;
    yy(65,1) = -0.25 * F(65,1);
    yy(65,65) = -0.25 * F(65,65);
    yy = ReconfromAmpPhase(yy,Phase);
    
    for i = 1:128
        for j = 1:128
            D = [xx(i,j) xy(i,j); xy(i,j) yy(i,j)];
            e = eig(D);
            eigenvalues(i,j,1) = e(1);
            eigenvalues(i,j,2) = e(2);  
        end
    end
    
    for j = 1:65
        for k = 1:128
            if k < 66
                Fx(j,k) = F(j,k)*(k-1)/128;
            else
                Fx(j,k) = F(j,k)*(k-129)/128;
            end
        end
    end
  
    Fx(1,1) = 0;
    Fx(1,65) = 0;
    Fx(65,1) = 0;
    Fx(65,65) = 0;
    
    Phase = Phase + pi/2;
    Phase(1,1) = 0;
    Phase(1,65) = 0;
    Phase(65,1) = 0;
    Phase(65,65) = 0;
            
    IA = ReconfromAmpPhase(Fx,Phase);
    for j = 1:65
        for k = 1:128
            Fy(j,k) = F(j,k)*(j-1)/128;
        end
    end
    Fy(1,1) = 0;
    Fy(1,65) = 0;
    Fy(65,1) = 0;
    Fy(65,65) = 0;
    IB = ReconfromAmpPhase(Fy,Phase);
    
    for i = 1:128
        for j = 1:128
            if eigenvalues(i,j,1) > 0 && eigenvalues(i,j,2) > 0
                image_e(i,j) = 0;
            else
                D = [xx(i,j) xy(i,j); xy(i,j) yy(i,j)];
                [V,D] = eig(D);
                if eigenvalues(i,j,1) > eigenvalues(i,j,2)
                    u = V(:,2);
                else
                    u = V(:,1);
                end
                image_e(i,j) = abs([IA(i,j) IB(i,j)]*u);
            end
        end
    end
   
    image_e(1,:) = 0;
    image_e(:,1) = 0;
    image_e(128,:) = 0;
    image_e(:,128) = 0;
    imshow(image_e)
    count = 0;
    for i = 1:128
        for j = 1:128
            if image_e(i,j) <= 5
                R(i,j) = 0;
            else
                if eigenvalues(i,j,1) > eigenvalues(i,j,2)
                    u = V(:,2);
                else
                    u = V(:,1);
                end
                theta = atan(u(2)/u(1));
                if theta < -3*pi/8 || theta > 3*pi/8
                    if (i == 1 || image_e(i,j) > image_e(i-1,j)) && (i == 128 || image_e(i,j) > image_e(i+1,j))
                        R(i,j) = 1;
                    else
                        R(i,j) = 0;
                    end
                elseif theta > pi/8 && theta < 3*pi/8
                    if (i == 1 || j == 128 || image_e(i,j) > image_e(i-1,j+1)) && (i == 128 || j == 1 || image_e(i,j) > image_e(i+1,j-1))
                        R(i,j) = 1;
                    else
                        R(i,j) = 0;
                    end
                elseif theta < pi/8 && theta > -pi/8
                    if (j == 1 || image_e(i,j) > image_e(i,j-1)) && (j == 128 || image_e(i,j) > image_e(i,j+1))
                        R(i,j) = 1;
                    else
                        R(i,j) = 0;
                    end
                else
                    if (i == 128 || j == 1 || image_e(i,j) > image_e(i+1,j-1)) && (i == 1 || j == 128 || image_e(i,j) > image_e(i-1,j+1))
                        R(i,j) = 1;
                    else
                        R(i,j) = 0;
                    end
                end
            end
        end
        
    end
    imshow(R)
end

function [out] = binomial_blurring(img)


    [h, w] = size(img); 
    out = im2double(img);
    for i = 1:h
        for j = 1:w
            if j == w
                out(i,j) = out(i,j)*2 / 2;
            else
                out(i,j) = (out(i,j) + out(i, j+1)) / 2;
            end
        end
        for j = w:-1:1
            if j == 1
                out(i,j) = out(i,j)*2 / 2;
            else
                out(i,j) = (out(i,j) + out(i, j-1)) / 2;
            end
        end
    end

        
    for j = 1:w
        for i = 1:h
            if i == h
                out(i,j) = out(i,j)*2 / 2;
            else
                out(i,j) = (out(i,j) + out(i+1, j)) / 2;
            end
        end

        for i = h:-1:1
            if i == 1
                out(i,j) = out(i,j)*2 / 2;
            else
                out(i,j) = (out(i,j) + out(i-1, j)) / 2;
            end
        end
    end
        
end
function noise_level = gnoise(sigma)
% try using 0.01 times the RMS intensity to be the sigma value 
  noise_level = ((rand + rand + rand + rand + rand + rand + rand + rand +rand + rand + rand + rand) - 6)*sigma; 
end
function [IA, std_dev, IB] = problemOne(img)
    A1 = img(:,:,1);
    A2 = img(:,:,2);
    A3 = img(:,:,3);
    A1 = im2double(A1);
    A2 = im2double(A2);
    A3 = im2double(A3);
    img1 = sqrt(A1.^2 + A2.^2 + A3.^2);
    [x1,y1] = size(img1);
    sep_x1 = floor(x1/128);
    sep_y1 = floor(y1/128);
    for i = 1:128
        for j = 1:128
            IA(i,j) = img1(i*sep_x1,j*sep_y1);
        end
    end
    for j = 1:128
        for k = 1:128
            if j < k
                IB(j,k) = 1000;
            end
            if j > k
                IB(j,k) = -1000;
            end
            if j == k
                IB(j,k) = 0;
            end
        end
    end
    IB = binomial_blurring(IB);
    sum = 0;
    for i = 1:128
        for j = 1:128
            sum = sum + IB(i,j)^2;
        end
    end
    std_dev = 0.01*sqrt(sum/128^2);
   
    for i = 1:128
        for j = 1:128
            noise_level = gnoise(std_dev);
            IB(i,j) = IB(i,j) + noise_level;
        end
    end
    

        
    % IA = 128x128 sub-image
    % std_dev = standard deviation (0.01 times the RMS intensity of the blurred diagonal image)
    % IB = 128x128 diagonal image with binomial_blurring and gnoise

end
    
