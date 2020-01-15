##################################################
# Polyacetylene
##################################################

# Find the actual crossing point
x0, y0 = TheoryOfCISS.find_crossing(
    Polyacetylene(N = 38),
    xsymbol = :θ,
    x0 = 1.5,
    y0 = 4.9,
    threshold=1e-6,
    attenuation=1e-2,
    correct_overlaps = false
);

# Calculate the D and C vectors
vs = TheoryOfCISS.atpoint(
    polarization,
    x0,
    y0,
    Polyacetylene(N = 38),
    xsymbol = :θ,
    ysymbol = :E,
    calc = TheoryOfCISS.calc_data1,
    correct_overlaps = false,
);

##################################################
# Helicene
##################################################

let
    f = polarization
    xsymbol = :δz
    x0 = 0.8
    ysymbol = :E
    y0 = -3.5
    mol = Helicene(N=7)

    x0, y0 = TheoryOfCISS.find_crossing(
        x0, y0, mol,
        xsymbol = xsymbol,
        threshold=1e-6,
        attenuation=1e-2
    )

    TheoryOfCISS.atpoint(
        f, x0, y0, mol,
        xsymbol = xsymbol,
        ysymbol = ysymbol,
        calc = TheoryOfCISS.calc_data1,
    )
end
