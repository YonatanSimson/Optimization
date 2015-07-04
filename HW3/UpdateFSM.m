function [distMap] =UpdateFSM(distMap,M11_1,M21_1,M12_1,M22_1,M11_2,M21_2,M12_2,M22_2,g11,g22,g12,m,n)


    %Direction 3 - top to bottom
        %FMM update

        for i=2:(m-1)
            G11=g11(i,2:end-1);
            G22=g22(i,2:end-1);
            G12=g12(i,2:end-1);
           %---------------------------------from triangle 1-------------------------------------
            t1_1=distMap(i+M21_1,1:end-2);
            t2_1=distMap(i+M22_1,2:end-1);

            E11_1=M11_1*(G11*M11_1 + G12*M21_1) + M21_1*(G12*M11_1 + G22*M21_1);
            E12_1=M12_1*(G11*M11_1 + G12*M21_1) + M22_1*(G12*M11_1 + G22*M21_1);
            E21_1=M11_1*(G11*M12_1 + G12*M22_1) + M21_1*(G12*M12_1 + G22*M22_1);
            E22_1=M12_1*(G11*M12_1 + G12*M22_1) + M22_1*(G12*M12_1 + G22*M22_1);
            

            determinant=(E11_1.*E22_1 - E12_1.*E21_1);
            Q11_1=E22_1./determinant;
            Q12_1=-E12_1./determinant;
            Q21_1=-E21_1./determinant;
            Q22_1=E11_1./determinant;


            tQt=t1_1.*(Q11_1.*t1_1 + Q12_1.*t2_1) + t2_1.*(Q12_1.*t1_1 + Q22_1.*t2_1);  %*********
            onesQones=Q11_1 + Q12_1 + Q21_1 + Q22_1;
            onesQt=t1_1.*(Q11_1 + Q21_1) + t2_1.*(Q12_1 + Q22_1);
            disc=(2*onesQt).^2-4*onesQones.*(tQt-1);
            %consistent=(disc>0).*(t1_1~=inf).*(t2_1~=inf);
            t01=(2*onesQt+ sqrt(abs(disc)))./(2*onesQones);
            t02=(2*onesQt- sqrt(abs(disc)))./(2*onesQones); 
            t0_newQuad=max(t01,t02); %t0_newQuad is the t0 found from the quadratic equation
            %checking monotonicity Q(t-t0)<0 or Q(t0-t)>0
            checkRow1= Q11_1.*(t0_newQuad - t1_1) + Q12_1.*(t0_newQuad - t2_1); %Q(t0-t) row 1
            checkRow2= Q21_1.*(t0_newQuad - t1_1) + Q22_1.*(t0_newQuad - t2_1); %Q(t0-t) row 2
                
            
            consistent=(disc>=0);
            isMonCon=(checkRow1>0).*(checkRow2>0).*(consistent);

            %Dijkstra update
            t0_newDijkstra1=t1_1+sqrt(M11_1*(G11*M11_1 + G12*M21_1) + M21_1*(G12*M11_1 + G22*M21_1));
            t0_newDijkstra2=t2_1+sqrt(M12_1*(G11*M12_1 + G12*M22_1) + M22_1*(G12*M12_1 + G22*M22_1));
            t0_newDijkstra=min( t0_newDijkstra1,t0_newDijkstra2);
            UpdatingNode=(t0_newDijkstra1>t0_newDijkstra2)+1;%      =1 if the updating node is 1 , =2 if the updating node is 2

            t0_new1= t0_newQuad.*isMonCon+t0_newDijkstra.*(1-isMonCon);
            DijkstraUpdate1=(1-isMonCon).*UpdatingNode; % =0 if the update was trough FM , =1 if the update is Dijkstra and the updating node is 1 , =2 if the update is Dijkstra and the updating node is 2


            %---------------------------------from triangle 2-------------------------------------

             %from triangle 2
            t1_2=distMap(i+M21_2,2:end-1);
            t2_2=distMap(i+M22_2,3:end);
            E11_2=M11_2*(G11*M11_2 + G12*M21_2) + M21_2*(G12*M11_2 + G22*M21_2);
            E12_2=M12_2*(G11*M11_2 + G12*M21_2) + M22_2*(G12*M11_2 + G22*M21_2);
            E21_2=M11_2*(G11*M12_2 + G12*M22_2) + M21_2*(G12*M12_2 + G22*M22_2);
            E22_2=M12_2*(G11*M12_2 + G12*M22_2) + M22_2*(G12*M12_2 + G22*M22_2);



            determinant=(E11_2.*E22_2 - E12_2.*E21_2);
            Q11_2=E22_2./determinant;
            Q12_2=-E12_2./determinant;
            Q21_2=-E21_2./determinant;
            Q22_2=E11_2./determinant;


            tQt=t1_2.*(Q11_2.*t1_2 + Q21_2.*t2_2) + t2_2.*(Q12_2.*t1_2 + Q22_2.*t2_2);
            onesQones=Q11_2 + Q12_2 + Q21_2 + Q22_2;
            onesQt=t1_2.*(Q11_2 + Q21_2) + t2_2.*(Q12_2 + Q22_2);
            disc=(2*onesQt).^2-4*onesQones.*(tQt-1);
            t01=(2*onesQt+ sqrt(abs(disc)))./(2*onesQones);
            t02=(2*onesQt- sqrt(abs(disc)))./(2*onesQones); 
            t0_newQuad=max(t01,t02); %t0_newQuad is the t0 found from the quadratic equation
            %checking monotonicity Q(t-t0)<0 or Q(t0-t)>0
            checkRow1= Q11_2.*(t0_newQuad - t1_2) + Q12_2.*(t0_newQuad - t2_2); %Q(t0-t) row 1
            checkRow2= Q21_2.*(t0_newQuad - t1_2) + Q22_2.*(t0_newQuad - t2_2); %Q(t0-t) row 2
                
            consistent=(disc>=0);
            isMonCon=(checkRow1>0).*(checkRow2>0).*(consistent);

            %Dijkstra update
            t0_newDijkstra1=t1_2+sqrt(M11_2*(G11*M11_2 + G12*M21_2) + M21_2*(G12*M11_2 + G22*M21_2));
            t0_newDijkstra2=t2_2+sqrt(M12_2*(G11*M12_2 + G12*M22_2) + M22_2*(G12*M12_2 + G22*M22_2));
            t0_newDijkstra=min( t0_newDijkstra1,t0_newDijkstra2);
            UpdatingNode=(t0_newDijkstra1>t0_newDijkstra2)+1;%      =1 if the updating node is 1 , =2 if the updating node is 2

            t0_new2= t0_newQuad.*isMonCon+t0_newDijkstra.*(1-isMonCon);
            DijkstraUpdate2=(1-isMonCon).*UpdatingNode; % =0 if the update was trough FM , =1 if the update is Dijkstra and the updating node is 1 , =2 if the update is Dijkstra and the updating node is 2




            %---------------------------------------------------------------Updating-------------------------------------------------
            t0_old=distMap(i,2:end-1);
            %t0_new=min(t0_new1,t0_new2);
            update1=(t0_new2-t0_new1>0);
            %update2=1-update1;
            t0_new=update1.*t0_new1+(1-update1).*t0_new2;
            
            update_new=(t0_old-t0_new)>0;
            distMap(i,2:end-1)=t0_new.*update_new+t0_old.*(1-update_new);
            
        end 
    
        
end