function y = polylog_BEC_fermi(z)
y = z + z.^2/(2*sqrt(2)) + z.^3/(3*sqrt(3)) + z.^4/8 + ...
z.^5/(5*sqrt(5)) + z.^6/(6*sqrt(6)) + z.^7/(7*sqrt(7)) +...
z.^8/(16*sqrt(2)) + z.^9/27 + z.^10/(10*sqrt(10));
end
