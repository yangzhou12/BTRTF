% BTRTF demo for image denoising
addpath(genpath('BTRTF'));

% 'Test_Images' consists of 8 sample images from the BSD500 dataset
imgDir = './Test_Images';
D= dir(fullfile(imgDir,'*.jpg'));

for i = 1:numel(D)   
    % Load images
    imgName = fullfile(imgDir,D(i).name);
    Img = double(imread(imgName));
    
    X_true = Img/255;
    Is = size(X_true);
    maxFea = max(X_true(:));
    
    % Generate sparse outliers
    SparseRatio = 0.10;
    E_true = zeros(size(X_true));    
    Eomega = randsample(prod(Is), round(prod(Is)*SparseRatio));
    E_true(Eomega) = rand(length(Eomega),1); % the sparse component

    % Observation
    Y = X_true + E_true;
    
    % BTRTF training
    model = BTRTF( Y, 'init', 'ml', 'Rs', 150*ones(Is(3),1), ...
        'initVar', 1e-7, 'maxIters', 500, 'tol', 1e-4 );
    Xhat = model.TX;
    
    imshow(Xhat); % Show the recovered image
    PSNR = myPSNR(X_true, Xhat, maxFea);
    fprintf('%d-th image, PSNR:%.4f \n', i, PSNR);
end
    