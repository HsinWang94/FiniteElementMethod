module SpaceTruss

    export SpaceTrussAssemble,SpaceTrussElementForce,SpaceTrussElementLength,SpaceTrussElementStiffness,SpaceTrussElementStress

    function SpaceTrussAssemble(K,k,i,j)
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

    function SpaceTrussElementForce(E,A,L,thetax,thetay,thetaz,u)
        x = thetax * pi/180;
        w = thetay * pi/180;
        v = thetaz * pi/180;
        Cx = cos(x);
        Cy = cos(w);
        Cz = cos(v);
        y = E*A/L*[-Cx -Cy -Cz Cx Cy Cz]*u;
        return y
    end

    function SpaceTrussElementLength(x1,y1,z1,x2,y2,z2)
        y = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
        return y
    end

    function SpaceTrussElementStress(E,L,thetax,thetay,thetaz,u)
        x = thetax * pi/180;
        w = thetay * pi/180;
        v = thetaz * pi/180;
        Cx = cos(x);
        Cy = cos(w);
        Cz = cos(v);
        y = E/L*[-Cx -Cy -Cz Cx Cy Cz]*u;
        return y
    end

    function SpaceTrussElementStiffness(E,A,L,thetax,thetay,thetaz)
        x = thetax*pi/180;
        u = thetay*pi/180;
        v = thetaz*pi/180;
        Cx = cos(x);
        Cy = cos(u);
        Cz = cos(v);
        w = [Cx*Cx Cx*Cy Cx*Cz ; Cy*Cx Cy*Cy Cy*Cz ; Cz*Cx Cz*Cy Cz*Cz];
        y = E*A/L*[w -w ; -w w];
        return y
    end
end