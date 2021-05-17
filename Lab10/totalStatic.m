function ptotal = totalStatic(M, gamma, pstatic)


ptotal = pstatic*((1 + ((gamma - 1)*0.5)*(M^2))^(gamma/(gamma - 1)));


end