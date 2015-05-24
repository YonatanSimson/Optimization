function IPlus = ActiveSet(x, gradf_x, lb, ub, epsilon)

IPlus = (gradf_x > 0 & lb <= x           & x <= lb + epsilon) | ...
        (gradf_x < 0 & ub - epsilon <= x & x <= ub);
