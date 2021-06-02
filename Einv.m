function xprime = Einv(x)
    global E
    xprime = floor(interp1(E, (1:length(E)), x));
end