function IPlus = ActiveSet(gradf_x, epsilon)

IPlus = abs(gradf_x) > epsilon;
