%%%%%%%%%%%%%%%%%%%%%%%--------------------------------
%
%
%------------------------------------------------------------------

function  Image = Zeropadding(input,Width,Hight)

S = size(input);
if( S(1)>Width||S(2)>Hight)

    error('Wrong use of the function!!!!!!!');
    
end

Width = round(Width);
Hight = round(Hight);

Image = zeros(Width,Hight);
Image(round(Width/2 - S(1)/2) + 1:round(Width/2 + S(1)/2) , round(Hight/2 - S(2)/2 + 1):round(Hight/2 + S(2)/2) ) = input;

end