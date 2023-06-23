#LinearBar.jl
module LinearBar

    export LinearBarAssemble,LinearBarElementForces,LinearBarElementStiffness,LinearBarElementStresses
"""

"""
    function LinearBarElementStiffness(E,A,L)    
        y = [E*A/L -E*A/L ; -E*A/L E*A/L];
        return y
    end
"""

"""
    function LinearBarAssemble(K,k,i,j)
        K[i,i] = K[i,i] + k[1,1];
        K[i,j] = K[i,j] + k[1,2];
        K[j,i] = K[j,i] + k[2,1];
        K[j,j] = K[j,j] + k[2,2];
        return K
    end
"""

"""
    function LinearBarElementForces(k,u)
        return k * u;
    end
"""

"""
    function LinearBarElementStresses(k, u, A)
        return k * u/A;
    end

end #module