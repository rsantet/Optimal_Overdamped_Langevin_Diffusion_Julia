"""Examples of potential energy functions defined on the torus.
"""

function cos_2(x)
    return cos(2 * π * x)
end

function cos_4(x)
    return cos(4 * π * x)
end

function cos_8(x)
    return cos(8 * π * x)
end

function cos_16(x)
    return cos(16 * π * x)
end

function sin_squared_cubed(x)
    return (
        sin(π / 2 + π * x)^2
        +
        sin(π / 2 + 2 * π * (x - 0.2))^3
    )
end

function sin_two_wells(x)
    return sin(4 * π * x) * (2 + sin(2 * π * x))
end

function sin_four_wells(x)
    return sin(8 * π * x) * (2 + sin(4 * π * x))
end