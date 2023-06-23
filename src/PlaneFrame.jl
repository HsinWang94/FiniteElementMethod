module PlaneFrame

    export PlaneFrameAssemble,PlaneFrameElementForces,PlaneFrameElementLength,PlaneFrameElementStiffness,PlaneFrameInclinedSupport

    function PlaneFrameAssemble(K,k,i,j)
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

    function PlaneFrameElementForces(E,A,I,L,theta,u)()
        x = theta * pi/180;
        C = cos(x);
        S = sin(x);
        w1 = E*A/L;
        w2 = 12*E*I/(L*L*L);
        w3 = 6*E*I/(L*L);
        w4 = 4*E*I/L;
        w5 = 2*E*I/L;
        kprime = [w1 0 0 -w1 0 0 ; 0 w2 w3 0 -w2 w3 ;
           0 w3 w4 0 -w3 w5 ; -w1 0 0 w1 0 0 ;
           0 -w2 -w3 0 w2 -w3 ; 0 w3 w5 0 -w3 w4];
        T = [C S 0 0 0 0 ; -S C 0 0 0 0 ; 0 0 1 0 0 0 ;
           0 0 0 C S 0 ; 0 0 0 -S C 0 ; 0 0 0 0 0 1];
        y = kprime*T* u;
        return y        
    end

    function PlaneFrameElementLength(x1,y1,x2,y2)
        y = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        return y        
    end

    function PlaneFrameElementStiffness(E,A,I,L,theta)
        x = theta*pi/180;
        C = cos(x);
        S = sin(x);
        w1 = A*C*C + 12*I*S*S/(L*L);
        w2 = A*S*S + 12*I*C*C/(L*L);
        w3 = (A-12*I/(L*L))*C*S;
        w4 = 6*I*S/L;
        w5 = 6*I*C/L;
        y = E/L*[w1 w3 -w4 -w1 -w3 -w4 ; w3 w2 w5 -w3 -w2 w5 ;
        -w4 w5 4*I w4 -w5 2*I ; -w1 -w3 w4 w1 w3 w4 ;
        -w3 -w2 -w5 w3 w2 -w5 ; -w4 w5 2*I w4 -w5 4*I];
        return y
    end

    function PlaneFrameInclinedSupport(T,i,alpha)
        x = alpha*pi/180;
        T[3*i-2,3*i-2] = cos(x);
        T[3*i-2,3*i-1] = sin(x);
        T[3*i-2,3*i] = 0;
        T[3*i-1,3*i-2] = -sin(x);
        T[3*i-1,3*i-1] = cos(x);
        T[3*i-1,3*i] = 0;
        T[3*i,3*i-2] = 0;
        T[3*i,3*i-1] = 0;
        T[3*i,3*i] = 1;
        y = T;
        return
    end
end