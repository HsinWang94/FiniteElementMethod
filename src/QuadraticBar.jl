#QuadraticBar.jl
module QuadraticBar

    export QuadraticBarElementStiffness,QuadraticBarAssemble,QuadraticBarElementForces,QuadraticBarElementStresses
"""
"""
    function QuadraticBarElementStiffness(E,A,L)
        y = E*A/(3*L)*[7 1 -8 ; 1 7 -8 ; -8 -8 16];
        return y
    end
"""
"""
    function QuadraticBarAssemble(K,k,i,j,m)
        K[i,i] = K[i,i] + k[1,1];
        K[i,j] = K[i,j] + k[1,2];
        K[i,m] = K[i,m] + k[1,3];
        K[j,i] = K[j,i] + k[2,1];
        K[j,j] = K[j,j] + k[2,2];
        K[j,m] = K[j,m] + k[2,3];
        K[m,i] = K[m,i] + k[3,1];
        K[m,j] = K[m,j] + k[3,2];
        K[m,m] = K[m,m] + k[3,3];
        return K
    end
"""
"""
    function QuadraticBarElementForces(k,u)
        y = k * u;
        return y
    end
"""
"""
    function QuadraticBarElementStresses(k, u, A)
        y = k * u/A;
        return y
    end
    
end