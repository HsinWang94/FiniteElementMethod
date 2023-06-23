#Spring.jl
module Spring

    export SpringElementStiffness,SpringAssemble,SpringElementForces

"""

"""
function SpringElementStiffness(k)
    y = [k -k ; -k k];
    return y
end
"""

"""
function SpringAssemble(K,k,i,j)
    K[i,i] = K[i,i] + k[1,1];
    K[i,j] = K[i,j] + k[1,2];
    K[j,i] = K[j,i] + k[2,1];
    K[j,j] = K[j,j] + k[2,2];
    return K
end
"""

"""
function SpringElementForces(k,u)
    return k * u;
end

end #module