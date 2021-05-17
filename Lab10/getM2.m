function M2 = getM2(M1, gamma)

top = 1 + (gamma - 1)*0.5*(M1^2);
bot = gamma*(M1^2) - (gamma - 1)*0.5;

M2 = sqrt(top/bot);

end