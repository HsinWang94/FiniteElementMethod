module Grid

    export GridAssemble,GridElementForces,GridElementLength,GridElementStiffness

    function GridAssemble(K,k,i,j)
        K[3*i-2,3*i-2] = K[3*i-2,3*i-2] + k[1,1];
        K[3*i-2,3*i-1] = K[3*i-2,3*i-1] + k[1,2];
        K[3*i-2,3*i] = K[3*i-2,3*i] + k[1,3];
        K[3*i-2,3*j-2] = K[3*i-2,3*j-2] + k[1,4];
        K[3*i-2,3*j-1] = K[3*i-2,3*j-1] + k[1,5];
        K[3*i-2,3*j] = K[3*i-2,3*j] + k[1,6];
        K[3*i-1,3*i-2] = K[3*i-1,3*i-2] + k[2,1];
        K[3*i-1,3*i-1] = K[3*i-1,3*i-1] + k[2,2];
        K[3*i-1,3*i] = K[3*i-1,3*i] + k[2,3];
        K[3*i-1,3*j-2] = K[3*i-1,3*j-2] + k[2,4];
        K[3*i-1,3*j-1] = K[3*i-1,3*j-1] + k[2,5];
        K[3*i-1,3*j] = K[3*i-1,3*j] + k[2,6];
        K[3*i,3*i-2] = K[3*i,3*i-2] + k[3,1];
        K[3*i,3*i-1] = K[3*i,3*i-1] + k[3,2];
        K[3*i,3*i] = K[3*i,3*i] + k[3,3];
        K[3*i,3*j-2] = K[3*i,3*j-2] + k[3,4];
        K[3*i,3*j-1] = K[3*i,3*j-1] + k[3,5];
        K[3*i,3*j] = K[3*i,3*j] + k[3,6];
        K[3*j-2,3*i-2] = K[3*j-2,3*i-2] + k[4,1];
        K[3*j-2,3*i-1] = K[3*j-2,3*i-1] + k[4,2];
        K[3*j-2,3*i] = K[3*j-2,3*i] + k[4,3];
        K[3*j-2,3*j-2] = K[3*j-2,3*j-2] + k[4,4];
        K[3*j-2,3*j-1] = K[3*j-2,3*j-1] + k[4,5];
        K[3*j-2,3*j] = K[3*j-2,3*j] + k[4,6];
        K[3*j-1,3*i-2] = K[3*j-1,3*i-2] + k[5,1];
        K[3*j-1,3*i-1] = K[3*j-1,3*i-1] + k[5,2];
        K[3*j-1,3*i] = K[3*j-1,3*i] + k[5,3];
        K[3*j-1,3*j-2] = K[3*j-1,3*j-2] + k[5,4];
        K[3*j-1,3*j-1] = K[3*j-1,3*j-1] + k[5,5];
        K[3*j-1,3*j] = K[3*j-1,3*j] + k[5,6];
        K[3*j,3*i-2] = K[3*j,3*i-2] + k[6,1];
        K[3*j,3*i-1] = K[3*j,3*i-1] + k[6,2];
        K[3*j,3*i] = K[3*j,3*i] + k[6,3];
        K[3*j,3*j-2] = K[3*j,3*j-2] + k[6,4];
        K[3*j,3*j-1] = K[3*j,3*j-1] + k[6,5];
        K[3*j,3*j] = K[3*j,3*j] + k[6,6];
        return K
    end

    function GridElementForces(E,G,I,J,L,theta,u)
        x = theta*pi/180;
        C = cos(x);
        S = sin(x);
        w1 = 12*E*I/(L*L*L);
        w2 = 6*E*I/(L*L);
        w3 = G*J/L;
        w4 = 4*E*I/L;
        w5 = 2*E*I/L;
        kprime = [w1 0 w2 -w1 0 w2 ; 0 w3 0 0 -w3 0 ;
        w2 0 w4 -w2 0 w5 ; -w1 0 -w2 w1 0 -w2 ;
        0 -w3 0 0 w3 0 ; w2 0 w5 -w2 0 w4];
        R = [1 0 0 0 0 0 ; 0 C S 0 0 0 ; 0 -S C 0 0 0 ;
        0 0 0 1 0 0 ; 0 0 0 0 C S ; 0 0 0 0 -S C];
        y = kprime*R* u;
        return y
    end

    function GridElementLength(x1,y1,x2,y2)
        y = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        return y        
    end

    function GridElementStiffness(E,G,I,J,L,theta)
        x = theta*pi/180;
        C = cos(x);
        S = sin(x);
        w1 = 12*E*I/(L*L*L);
        w2 = 6*E*I/(L*L);
        w3 = G*J/L;
        w4 = 4*E*I/L;
        w5 = 2*E*I/L;
        kprime = [w1 0 w2 -w1 0 w2 ; 0 w3 0 0 -w3 0 ;
           w2 0 w4 -w2 0 w5 ; -w1 0 -w2 w1 0 -w2 ;
           0 -w3 0 0 w3 0 ; w2 0 w5 -w2 0 w4];
        R = [1 0 0 0 0 0 ; 0 C S 0 0 0 ; 0 -S C 0 0 0 ;
           0 0 0 1 0 0 ; 0 0 0 0 C S ; 0 0 0 0 -S C];
        y = R'*kprime*R;
        return y        
    end
end