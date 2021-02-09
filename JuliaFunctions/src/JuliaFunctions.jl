module JuliaFunctions
    using BlackBoxOptim
    function f(x,y)
        z = 0
        for i in 1:x
            z += i ^ y
        end
        return z
    end
    export f
end  # module
