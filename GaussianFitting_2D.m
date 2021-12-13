
%% ---------User Input---------------------
function [Output FitPara]= GaussianFitting_2D( Input , Initial , FigureID )
  
   
    MdataSize = size(Input , 1); % Size of nxn data matrix
    % parameters are: [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]
    x0 = Initial; %Inital guess parameters

    %% ---Generate centroid to be fitted--------------------------------------
    HSize = floor( MdataSize/2 );
    [X,Y] = meshgrid(-HSize:HSize);
    xdata = zeros(size(X,1),size(Y,2),2);
    xdata(:,:,1) = X;
    xdata(:,:,2) = Y;
    [Xhr,Yhr] = meshgrid(linspace(-HSize,HSize,HSize*2+1)); % generate high res grid for plot
    xdatahr = zeros( HSize*2+1 , HSize*2+1 , 2 );
    xdatahr(:,:,1) = Xhr;
    xdatahr(:,:,2) = Yhr;

    %% --- Fit---------------------
        x0 =x0(1:5);
        xin(6) = 0; 
     
        lb = [0,-HSize,0,-HSize,0];
        ub = [realmax('double'),HSize,(HSize)^2,HSize,(HSize)^2];
        [FitPara,resnorm,residual,exitflag] =  lsqcurvefit(@D2GaussFunction,x0,xdata,Input,lb,ub);
  


    %% ---------Plot 3D Image-------------
    Output  = D2GaussFunction(FitPara,xdatahr);
%     figure( FigureID )
%     surface(Xhr,Yhr,Output,'EdgeColor','none') %plot fit
%        alpha(0.9)  

