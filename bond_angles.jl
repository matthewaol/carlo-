function calculate_bond_angles(coords)
    n_atoms = length(coords)
    for i in 2:(n_atoms - 1)
        angle = Molly.bond_angle(coords[i-1], coords[i], coords[i+1], boundary)
        println("Bond angle at atom $i: $angle radians")
    end
end