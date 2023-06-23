module Beam
    
    export BeamAssemble,BeamElementForces,BeamElementStiffness

    function BeamAssemble(K,k,i,j)
        K[2*i-1,2*i-1] = K[2*i-1,2*i-1] + k[1,1];
        K[2*i-1,2*i] = K[2*i-1,2*i] + k[1,2];
        K[2*i-1,2*j-1] = K[2*i-1,2*j-1] + k[1,3];
        K[2*i-1,2*j] = K[2*i-1,2*j] + k[1,4];
        K[2*i,2*i-1] = K[2*i,2*i-1] + k[2,1];
        K[2*i,2*i] = K[2*i,2*i] + k[2,2];
        K[2*i,2*j-1] = K[2*i,2*j-1] + k[2,3];
        K[2*i,2*j] = K[2*i,2*j] + k[2,4];
        K[2*j-1,2*i-1] = K[2*j-1,2*i-1] + k[3,1];
        K[2*j-1,2*i] = K[2*j-1,2*i] + k[3,2];
        K[2*j-1,2*j-1] = K[2*j-1,2*j-1] + k[3,3];
        K[2*j-1,2*j] = K[2*j-1,2*j] + k[3,4];
        K[2*j,2*i-1] = K[2*j,2*i-1] + k[4,1];
        K[2*j,2*i] = K[2*j,2*i] + k[4,2];
        K[2*j,2*j-1] = K[2*j,2*j-1] + k[4,3];
        K[2*j,2*j] = K[2*j,2*j] + k[4,4];
        return K
    end

    function BeamElementForces(k,u)
        y = k * u;
        return y       
    end

    function BeamElementStiffness(E,I,L)
        y = E*I/(L*L*L)*[12 6*L -12 6*L ; 6*L 4*L*L -6*L 2*L*L ;
        -12 -6*L 12 -6*L ; 6*L 2*L*L -6*L 4*L*L];
        return y
    end
end