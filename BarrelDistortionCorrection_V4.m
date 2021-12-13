function outuput = BarrelDistortionCorrection_V4( Input , K1 , K2, K3 , padding);



    Size_X = size(Input,1)+padding;
    Size_Y = size(Input,2)+padding;
    Input_Z = Zeropadding( Input , 2*Size_X+1 , 2*Size_Y+1 );


X         = -Size_X:1:Size_X;
Y         = -Size_Y:1:Size_Y;
[XX , YY] = meshgrid(X , Y);

XX_C      = XX.*(  1 + K1*XX.^2 + K2*YY.^2 + K3*(XX.^2+YY.^2).^2   );
YY_C      = YY.*(  1 + K1*XX.^2 + K2*YY.^2 + K3*(XX.^2+YY.^2).^2 );
XX_C      = round(XX_C + Size_X);
XX_C      = XX_C.*(XX_C>0).*(XX_C<=2*Size_X+1);

YY_C      = round(YY_C + Size_Y);
YY_C      = YY_C.*(YY_C>0).*(YY_C<=2*Size_Y+1);

temp = zeros(2*Size_X+1,2*Size_Y+1);
for Row = 1:2*Size_X+1
    for Col = 1:2*Size_Y+1

        if( XX_C(Row,Col)==0||YY_C(Row,Col)==0 )
            temp(Row , Col) = 0;
        else
            temp(Row , Col) = Input_Z( XX_C(Row,Col) , YY_C(Row,Col) );
        end
        
    end
end

outuput = temp( Size_X+1-floor(Size_X*0.5):Size_X+1+floor(Size_X*0.5) , Size_Y+1-floor(Size_Y*0.5):Size_Y+1+floor(Size_Y*0.5) );
outuput = transpose(outuput);

   
end


