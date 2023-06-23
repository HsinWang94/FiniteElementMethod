module SpaceFrame

    export SpaceFrameAssemble,SpaceFrameElementForces,SpaceFrameElementLength,SpaceFrameElementStiffness

    function SpaceFrameAssemble(K,k,i,j)
        K[6*i-5,6*i-5] = K[6*i-5,6*i-5] + k[1,1];
        K[6*i-5,6*i-4] = K[6*i-5,6*i-4] + k[1,2];
        K[6*i-5,6*i-3] = K[6*i-5,6*i-3] + k[1,3];
        K[6*i-5,6*i-2] = K[6*i-5,6*i-2] + k[1,4];
        K[6*i-5,6*i-1] = K[6*i-5,6*i-1] + k[1,5];
        K[6*i-5,6*i] = K[6*i-5,6*i] + k[1,6];
        K[6*i-5,6*j-5] = K[6*i-5,6*j-5] + k[1,7];
        K[6*i-5,6*j-4] = K[6*i-5,6*j-4] + k[1,8];
        K[6*i-5,6*j-3] = K[6*i-5,6*j-3] + k[1,9];
        K[6*i-5,6*j-2] = K[6*i-5,6*j-2] + k[1,10];
        K[6*i-5,6*j-1] = K[6*i-5,6*j-1] + k[1,11];
        K[6*i-5,6*j] = K[6*i-5,6*j] + k[1,12];
        K[6*i-4,6*i-5] = K[6*i-4,6*i-5] + k[2,1];
        K[6*i-4,6*i-4] = K[6*i-4,6*i-4] + k[2,2];
        K[6*i-4,6*i-3] = K[6*i-4,6*i-3] + k[2,3];
        K[6*i-4,6*i-2] = K[6*i-4,6*i-2] + k[2,4];
        K[6*i-4,6*i-1] = K[6*i-4,6*i-1] + k[2,5];
        K[6*i-4,6*i] = K[6*i-4,6*i] + k[2,6];
        K[6*i-4,6*j-5] = K[6*i-4,6*j-5] + k[2,7];
        K[6*i-4,6*j-4] = K[6*i-4,6*j-4] + k[2,8];
        K[6*i-4,6*j-3] = K[6*i-4,6*j-3] + k[2,9];
        K[6*i-4,6*j-2] = K[6*i-4,6*j-2] + k[2,10];
        K[6*i-4,6*j-1] = K[6*i-4,6*j-1] + k[2,11];
        K[6*i-4,6*j] = K[6*i-4,6*j] + k[2,12];
        K[6*i-3,6*i-5] = K[6*i-3,6*i-5] + k[3,1];
        K[6*i-3,6*i-4] = K[6*i-3,6*i-4] + k[3,2];
        K[6*i-3,6*i-3] = K[6*i-3,6*i-3] + k[3,3];
        K[6*i-3,6*i-2] = K[6*i-3,6*i-2] + k[3,4];
        K[6*i-3,6*i-1] = K[6*i-3,6*i-1] + k[3,5];
        K[6*i-3,6*i] = K[6*i-3,6*i] + k[3,6];
        K[6*i-3,6*j-5] = K[6*i-3,6*j-5] + k[3,7];
        K[6*i-3,6*j-4] = K[6*i-3,6*j-4] + k[3,8];
        K[6*i-3,6*j-3] = K[6*i-3,6*j-3] + k[3,9];
        K[6*i-3,6*j-2] = K[6*i-3,6*j-2] + k[3,10];
        K[6*i-3,6*j-1] = K[6*i-3,6*j-1] + k[3,11];
        K[6*i-3,6*j] = K[6*i-3,6*j] + k[3,12];
        K[6*i-2,6*i-5] = K[6*i-2,6*i-5] + k[4,1];
        K[6*i-2,6*i-4] = K[6*i-2,6*i-4] + k[4,2];
        K[6*i-2,6*i-3] = K[6*i-2,6*i-3] + k[4,3];
        K[6*i-2,6*i-2] = K[6*i-2,6*i-2] + k[4,4];
        K[6*i-2,6*i-1] = K[6*i-2,6*i-1] + k[4,5];
        K[6*i-2,6*i] = K[6*i-2,6*i] + k[4,6];
        K[6*i-2,6*j-5] = K[6*i-2,6*j-5] + k[4,7];
        K[6*i-2,6*j-4] = K[6*i-2,6*j-4] + k[4,8];
        K[6*i-2,6*j-3] = K[6*i-2,6*j-3] + k[4,9];
        K[6*i-2,6*j-2] = K[6*i-2,6*j-2] + k[4,10];
        K[6*i-2,6*j-1] = K[6*i-2,6*j-1] + k[4,11];
        K[6*i-2,6*j] = K[6*i-2,6*j] + k[4,12];
        K[6*i-1,6*i-5] = K[6*i-1,6*i-5] + k[5,1];
        K[6*i-1,6*i-4] = K[6*i-1,6*i-4] + k[5,2];
        K[6*i-1,6*i-3] = K[6*i-1,6*i-3] + k[5,3];
        K[6*i-1,6*i-2] = K[6*i-1,6*i-2] + k[5,4];
        K[6*i-1,6*i-1] = K[6*i-1,6*i-1] + k[5,5];
        K[6*i-1,6*i] = K[6*i-1,6*i] + k[5,6];
        K[6*i-1,6*j-5] = K[6*i-1,6*j-5] + k[5,7];
        K[6*i-1,6*j-4] = K[6*i-1,6*j-4] + k[5,8];
        K[6*i-1,6*j-3] = K[6*i-1,6*j-3] + k[5,9];
        K[6*i-1,6*j-2] = K[6*i-1,6*j-2] + k[5,10];
        K[6*i-1,6*j-1] = K[6*i-1,6*j-1] + k[5,11];
        K[6*i-1,6*j] = K[6*i-1,6*j] + k[5,12];
        K[6*i,6*i-5] = K[6*i,6*i-5] + k[6,1];
        K[6*i,6*i-4] = K[6*i,6*i-4] + k[6,2];
        K[6*i,6*i-3] = K[6*i,6*i-3] + k[6,3];
        K[6*i,6*i-2] = K[6*i,6*i-2] + k[6,4];
        K[6*i,6*i-1] = K[6*i,6*i-1] + k[6,5];
        K[6*i,6*i] = K[6*i,6*i] + k[6,6];
        K[6*i,6*j-5] = K[6*i,6*j-5] + k[6,7];
        K[6*i,6*j-4] = K[6*i,6*j-4] + k[6,8];
        K[6*i,6*j-3] = K[6*i,6*j-3] + k[6,9];
        K[6*i,6*j-2] = K[6*i,6*j-2] + k[6,10];
        K[6*i,6*j-1] = K[6*i,6*j-1] + k[6,11];
        K[6*i,6*j] = K[6*i,6*j] + k[6,12];
        K[6*j-5,6*i-5] = K[6*j-5,6*i-5] + k[7,1];
        K[6*j-5,6*i-4] = K[6*j-5,6*i-4] + k[7,2];
        K[6*j-5,6*i-3] = K[6*j-5,6*i-3] + k[7,3];
        K[6*j-5,6*i-2] = K[6*j-5,6*i-2] + k[7,4];
        K[6*j-5,6*i-1] = K[6*j-5,6*i-1] + k[7,5];
        K[6*j-5,6*i] = K[6*j-5,6*i] + k[7,6];
        K[6*j-5,6*j-5] = K[6*j-5,6*j-5] + k[7,7];
        K[6*j-5,6*j-4] = K[6*j-5,6*j-4] + k[7,8];
        K[6*j-5,6*j-3] = K[6*j-5,6*j-3] + k[7,9];
        K[6*j-5,6*j-2] = K[6*j-5,6*j-2] + k[7,10];
        K[6*j-5,6*j-1] = K[6*j-5,6*j-1] + k[7,11];
        K[6*j-5,6*j] = K[6*j-5,6*j] + k[7,12];
        K[6*j-4,6*i-5] = K[6*j-4,6*i-5] + k[8,1];
        K[6*j-4,6*i-4] = K[6*j-4,6*i-4] + k[8,2];
        K[6*j-4,6*i-3] = K[6*j-4,6*i-3] + k[8,3];
        K[6*j-4,6*i-2] = K[6*j-4,6*i-2] + k[8,4];
        K[6*j-4,6*i-1] = K[6*j-4,6*i-1] + k[8,5];
        K[6*j-4,6*i] = K[6*j-4,6*i] + k[8,6];
        K[6*j-4,6*j-5] = K[6*j-4,6*j-5] + k[8,7];
        K[6*j-4,6*j-4] = K[6*j-4,6*j-4] + k[8,8];
        K[6*j-4,6*j-3] = K[6*j-4,6*j-3] + k[8,9];
        K[6*j-4,6*j-2] = K[6*j-4,6*j-2] + k[8,10];
        K[6*j-4,6*j-1] = K[6*j-4,6*j-1] + k[8,11];
        K[6*j-4,6*j] = K[6*j-4,6*j] + k[8,12];
        K[6*j-3,6*i-5] = K[6*j-3,6*i-5] + k[9,1];
        K[6*j-3,6*i-4] = K[6*j-3,6*i-4] + k[9,2];
        K[6*j-3,6*i-3] = K[6*j-3,6*i-3] + k[9,3];
        K[6*j-3,6*i-2] = K[6*j-3,6*i-2] + k[9,4];
        K[6*j-3,6*i-1] = K[6*j-3,6*i-1] + k[9,5];
        K[6*j-3,6*i] = K[6*j-3,6*i] + k[9,6];
        K[6*j-3,6*j-5] = K[6*j-3,6*j-5] + k[9,7];
        K[6*j-3,6*j-4] = K[6*j-3,6*j-4] + k[9,8];
        K[6*j-3,6*j-3] = K[6*j-3,6*j-3] + k[9,9];
        K[6*j-3,6*j-2] = K[6*j-3,6*j-2] + k[9,10];
        K[6*j-3,6*j-1] = K[6*j-3,6*j-1] + k[9,11];
        K[6*j-3,6*j] = K[6*j-3,6*j] + k[9,12];
        K[6*j-2,6*i-5] = K[6*j-2,6*i-5] + k[10,1];
        K[6*j-2,6*i-4] = K[6*j-2,6*i-4] + k[10,2];
        K[6*j-2,6*i-3] = K[6*j-2,6*i-3] + k[10,3];
        K[6*j-2,6*i-2] = K[6*j-2,6*i-2] + k[10,4];
        K[6*j-2,6*i-1] = K[6*j-2,6*i-1] + k[10,5];
        K[6*j-2,6*i] = K[6*j-2,6*i] + k[10,6];
        K[6*j-2,6*j-5] = K[6*j-2,6*j-5] + k[10,7];
        K[6*j-2,6*j-4] = K[6*j-2,6*j-4] + k[10,8];
        K[6*j-2,6*j-3] = K[6*j-2,6*j-3] + k[10,9];
        K[6*j-2,6*j-2] = K[6*j-2,6*j-2] + k[10,10];
        K[6*j-2,6*j-1] = K[6*j-2,6*j-1] + k[10,11];
        K[6*j-2,6*j] = K[6*j-2,6*j] + k[10,12];
        K[6*j-1,6*i-5] = K[6*j-1,6*i-5] + k[11,1];
        K[6*j-1,6*i-4] = K[6*j-1,6*i-4] + k[11,2];
        K[6*j-1,6*i-3] = K[6*j-1,6*i-3] + k[11,3];
        K[6*j-1,6*i-2] = K[6*j-1,6*i-2] + k[11,4];
        K[6*j-1,6*i-1] = K[6*j-1,6*i-1] + k[11,5];
        K[6*j-1,6*i] = K[6*j-1,6*i] + k[11,6];
        K[6*j-1,6*j-5] = K[6*j-1,6*j-5] + k[11,7];
        K[6*j-1,6*j-4] = K[6*j-1,6*j-4] + k[11,8];
        K[6*j-1,6*j-3] = K[6*j-1,6*j-3] + k[11,9];
        K[6*j-1,6*j-2] = K[6*j-1,6*j-2] + k[11,10];
        K[6*j-1,6*j-1] = K[6*j-1,6*j-1] + k[11,11];
        K[6*j-1,6*j] = K[6*j-1,6*j] + k[11,12];
        K[6*j,6*i-5] = K[6*j,6*i-5] + k[12,1];
        K[6*j,6*i-4] = K[6*j,6*i-4] + k[12,2];
        K[6*j,6*i-3] = K[6*j,6*i-3] + k[12,3];
        K[6*j,6*i-2] = K[6*j,6*i-2] + k[12,4];
        K[6*j,6*i-1] = K[6*j,6*i-1] + k[12,5];
        K[6*j,6*i] = K[6*j,6*i] + k[12,6];
        K[6*j,6*j-5] = K[6*j,6*j-5] + k[12,7];
        K[6*j,6*j-4] = K[6*j,6*j-4] + k[12,8];
        K[6*j,6*j-3] = K[6*j,6*j-3] + k[12,9];
        K[6*j,6*j-2] = K[6*j,6*j-2] + k[12,10];
        K[6*j,6*j-1] = K[6*j,6*j-1] + k[12,11];
        K[6*j,6*j] = K[6*j,6*j] + k[12,12];
        return K
    end

    function SpaceFrameElementForces(E,G,A,Iy,Iz,J,x1,y1,z1,x2,y2,z2,u)
        L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
        w1 = E*A/L;
        w2 = 12*E*Iz/(L*L*L);
        w3 = 6*E*Iz/(L*L);
        w4 = 4*E*Iz/L;
        w5 = 2*E*Iz/L;
        w6 = 12*E*Iy/(L*L*L);
        w7 = 6*E*Iy/(L*L);
        w8 = 4*E*Iy/L;
        w9 = 2*E*Iy/L;
        w10 = G*J/L;
        kprime = [w1 0 0 0 0 0 -w1 0 0 0 0 0 ;
           0 w2 0 0 0 w3 0 -w2 0 0 0 w3 ;
           0 0 w6 0 -w7 0 0 0 -w6 0 -w7 0 ;
           0 0 0 w10 0 0 0 0 0 -w10 0 0 ;
           0 0 -w7 0 w8 0 0 0 w7 0 w9 0 ;
           0 w3 0 0 0 w4 0 -w3 0 0 0 w5 ;
           -w1 0 0 0 0 0 w1 0 0 0 0 0 ;
           0 -w2 0 0 0 -w3 0 w2 0 0 0 -w3 ;
           0 0 -w6 0 w7 0 0 0 w6 0 w7 0 ;
           0 0 0 -w10 0 0 0 0 0 w10 0 0 ;
           0 0 -w7 0 w9 0 0 0 w7 0 w8 0 ;
           0 w3 0 0 0 w5 0 -w3 0 0 0 w4];
        if x1 == x2 & y1 == y2
           if z2 > z1
              Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
           else
              Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
           end
        else
           CXx = (x2-x1)/L;
            CYx = (y2-y1)/L;
            CZx = (z2-z1)/L;
            D = sqrt(CXx*CXx + CYx*CYx);
            CXy = -CYx/D;
            CYy = CXx/D;
            CZy = 0;
            CXz = -CXx*CZx/D;
            CYz = -CYx*CZx/D;
            CZz = D;
            Lambda = [CXx CYx CZx ; CXy CYy CZy ; CXz CYz CZz];
        end
        R = [Lambda zeros(3) zeros(3) zeros(3) ; 
           zeros(3) Lambda zeros(3) zeros(3) ;
           zeros(3) zeros(3) Lambda zeros(3) ;
           zeros(3) zeros(3) zeros(3) Lambda];
        y = kprime*R* u;
        return y
    end

    function SpaceFrameElementLength(x1,y1,z1,x2,y2,z2)
        y = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
        return y
    end

    function SpaceFrameElementStiffness(E,G,A,Iy,Iz,J,x1,y1,z1,x2,y2,z2)
        L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
        w1 = E*A/L;
        w2 = 12*E*Iz/(L*L*L);
        w3 = 6*E*Iz/(L*L);
        w4 = 4*E*Iz/L;
        w5 = 2*E*Iz/L;
        w6 = 12*E*Iy/(L*L*L);
        w7 = 6*E*Iy/(L*L);
        w8 = 4*E*Iy/L;
        w9 = 2*E*Iy/L;
        w10 = G*J/L;
        kprime = [w1 0 0 0 0 0 -w1 0 0 0 0 0 ;
           0 w2 0 0 0 w3 0 -w2 0 0 0 w3 ;
           0 0 w6 0 -w7 0 0 0 -w6 0 -w7 0 ;
           0 0 0 w10 0 0 0 0 0 -w10 0 0 ;
           0 0 -w7 0 w8 0 0 0 w7 0 w9 0 ;
           0 w3 0 0 0 w4 0 -w3 0 0 0 w5 ;
           -w1 0 0 0 0 0 w1 0 0 0 0 0 ;
           0 -w2 0 0 0 -w3 0 w2 0 0 0 -w3 ;
           0 0 -w6 0 w7 0 0 0 w6 0 w7 0 ;
           0 0 0 -w10 0 0 0 0 0 w10 0 0 ;
           0 0 -w7 0 w9 0 0 0 w7 0 w8 0 ;
           0 w3 0 0 0 w5 0 -w3 0 0 0 w4];
        if x1 == x2 & y1 == y2
           if z2 > z1
              Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
           else
              Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
           end
        else
           CXx = (x2-x1)/L;
            CYx = (y2-y1)/L;
            CZx = (z2-z1)/L;
            D = sqrt(CXx*CXx + CYx*CYx);
            CXy = -CYx/D;
            CYy = CXx/D;
            CZy = 0;
            CXz = -CXx*CZx/D;
            CYz = -CYx*CZx/D;
            CZz = D;
            Lambda = [CXx CYx CZx ; CXy CYy CZy ; CXz CYz CZz];
        end
        R = [Lambda zeros(3) zeros(3) zeros(3) ; 
           zeros(3) Lambda zeros(3) zeros(3) ;
           zeros(3) zeros(3) Lambda zeros(3) ;
           zeros(3) zeros(3) zeros(3) Lambda];
        y = R'*kprime*R;
        return y        
    end
end