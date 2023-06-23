module LinearTriangle

    export LinearTriangleAssemble,LinearTriangleElementArea,LinearTriangleElementPStresses,LinearTriangleElementStiffness,LinearTriangleElementStresses
    
    function LinearTriangleAssemble(K,k,i,j,m)
        K[2*i-1,2*i-1] = K[2*i-1,2*i-1] + k[1,1];
        K[2*i-1,2*i] = K[2*i-1,2*i] + k[1,2];
        K[2*i-1,2*j-1] = K[2*i-1,2*j-1] + k[1,3];
        K[2*i-1,2*j] = K[2*i-1,2*j] + k[1,4];
        K[2*i-1,2*m-1] = K[2*i-1,2*m-1] + k[1,5];
        K[2*i-1,2*m] = K[2*i-1,2*m] + k[1,6];
        K[2*i,2*i-1] = K[2*i,2*i-1] + k[2,1];
        K[2*i,2*i] = K[2*i,2*i] + k[2,2];
        K[2*i,2*j-1] = K[2*i,2*j-1] + k[2,3];
        K[2*i,2*j] = K[2*i,2*j] + k[2,4];
        K[2*i,2*m-1] = K[2*i,2*m-1] + k[2,5];
        K[2*i,2*m] = K[2*i,2*m] + k[2,6];
        K[2*j-1,2*i-1] = K[2*j-1,2*i-1] + k[3,1];
        K[2*j-1,2*i] = K[2*j-1,2*i] + k[3,2];
        K[2*j-1,2*j-1] = K[2*j-1,2*j-1] + k[3,3];
        K[2*j-1,2*j] = K[2*j-1,2*j] + k[3,4];
        K[2*j-1,2*m-1] = K[2*j-1,2*m-1] + k[3,5];
        K[2*j-1,2*m] = K[2*j-1,2*m] + k[3,6];
        K[2*j,2*i-1] = K[2*j,2*i-1] + k[4,1];
        K[2*j,2*i] = K[2*j,2*i] + k[4,2];
        K[2*j,2*j-1] = K[2*j,2*j-1] + k[4,3];
        K[2*j,2*j] = K[2*j,2*j] + k[4,4];
        K[2*j,2*m-1] = K[2*j,2*m-1] + k[4,5];
        K[2*j,2*m] = K[2*j,2*m] + k[4,6];
        K[2*m-1,2*i-1] = K[2*m-1,2*i-1] + k[5,1];
        K[2*m-1,2*i] = K[2*m-1,2*i] + k[5,2];
        K[2*m-1,2*j-1] = K[2*m-1,2*j-1] + k[5,3];
        K[2*m-1,2*j] = K[2*m-1,2*j] + k[5,4];
        K[2*m-1,2*m-1] = K[2*m-1,2*m-1] + k[5,5];
        K[2*m-1,2*m] = K[2*m-1,2*m] + k[5,6];
        K[2*m,2*i-1] = K[2*m,2*i-1] + k[6,1];
        K[2*m,2*i] = K[2*m,2*i] + k[6,2];
        K[2*m,2*j-1] = K[2*m,2*j-1] + k[6,3];
        K[2*m,2*j] = K[2*m,2*j] + k[6,4];
        K[2*m,2*m-1] = K[2*m,2*m-1] + k[6,5];
        K[2*m,2*m] = K[2*m,2*m] + k[6,6];
        return K 
    end

    function LinearTriangleElementArea(xi,yi,xj,yj,xm,ym)
        y = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2;
        return y
    end

    function LinearTriangleElementPStresses(sigma)
        R = (sigma(1) + sigma(2))/2;
        Q = ((sigma(1) - sigma(2))/2)^2 + sigma(3)*sigma(3);
        M = 2*sigma(3)/(sigma(1) - sigma(2));
        s1 = R + sqrt(Q);
        s2 = R - sqrt(Q);
        theta = (atan(M)/2)*180/pi;
        y = [s1 ; s2 ; theta];
        return y
    end

    function LinearTriangleElementStiffness(E,NU,t,xi,yi,xj,yj,xm,ym,p)
        A = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2;
        betai = yj-ym;
        betaj = ym-yi;
        betam = yi-yj;
        gammai = xm-xj;
        gammaj = xi-xm;
        gammam = xj-xi;
        B = [betai 0 betaj 0 betam 0 ; 
        0 gammai 0 gammaj 0 gammam ;
        gammai betai gammaj betaj gammam betam]/(2*A);
        if p == 1 
            D = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2];
        elseif p == 2
            D = (E/(1+NU)/(1-2*NU))*[1-NU NU 0 ; NU 1-NU 0 ; 0 0 (1-2*NU)/2];
        end
        y = t*A*B'*D*B;
        return y
    end

    function LinearTriangleElementStresses(E,NU,t,xi,yi,xj,yj,xm,ym,p,u)
        A = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2;
        betai = yj-ym;
        betaj = ym-yi;
        betam = yi-yj;
        gammai = xm-xj;
        gammaj = xi-xm;
        gammam = xj-xi;
        B = [betai 0 betaj 0 betam 0 ; 
           0 gammai 0 gammaj 0 gammam ;
           gammai betai gammaj betaj gammam betam]/(2*A);
        if p == 1 
           D = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2];
        elseif p == 2
           D = (E/(1+NU)/(1-2*NU))*[1-NU NU 0 ; NU 1-NU 0 ; 0 0 (1-2*NU)/2];
        end
        y = D*B*u;
        return y
    end
end