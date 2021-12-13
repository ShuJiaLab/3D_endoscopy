function Output = Generate_RoundBeam( Radius , Image_Size )

         Output = zeros( Image_Size , Image_Size );
         Center = Image_Size/2;
         for rr = 1:Image_Size
             for cc=1:Image_Size
                 
                 
                 Dis = sqrt(  (rr-Center)^2 + (cc-Center)^2  );
                 
                 if( Dis<=Radius )
                     
                     Output(rr,cc)=255;
                 
             end
         end


end