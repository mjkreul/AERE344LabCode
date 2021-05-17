function Mma = getMma(ptotal, pstatic, gamma)

Mma = sqrt((((ptotal/pstatic)^((gamma - 1)/gamma)) - 1)*2/(gamma -1));

end