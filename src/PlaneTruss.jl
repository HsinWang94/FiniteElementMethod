# PlaneTruss.jl
module PlaneTruss

    export PlaneTrussAssemble,PlaneTrussElementLength,PlaneTrussElementForce,PlaneTrussElementStiffness,PlaneTrussElementStress,PlaneTrussInclinedSupport
    
    function PlaneTrussAssemble(K,k,i,j)
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

    function PlaneTrussElementLength(x1,y1,x2,y2)
        y = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        return y
    end
    
    function PlaneTrussElementForce(E,A,L,theta,u)
        x = theta * pi/180;
        C = cos(x);
        S = sin(x);
        y = E*A/L*[-C -S C S]* u;
        return y
    end
    
    function PlaneTrussElementStiffness(E,A,L, theta)
        x = theta*pi/180;
        C = cos(x);
        S = sin(x);
        y = E*A/L*[C*C C*S -C*C -C*S ; C*S S*S -C*S -S*S ;
           -C*C -C*S C*C C*S ; -C*S -S*S C*S S*S];
        return y
    end

    function PlaneTrussElementStress(E,L,theta,u)
        x = theta * pi/180;
        C = cos(x);
        S = sin(x);
        y = E/L*[-C -S C S]* u;
        return y
    end

    function PlaneTrussInclinedSupport(T,i,alpha)
        x = alpha*pi/180;
        T[2*i-1,2*i-1] = cos(x);
        T[2*i-1,2*i] = sin(x);
        T[2*i,2*i-1] = -sin(x) ;
        T[2*i,2*i] = cos(x);
        y = T;
        return y
    end    
end