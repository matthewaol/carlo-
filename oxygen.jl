using Molly
using GLMakie
using Statistics

function calculate_bond_angles(coords)
    bond_angles = []
    n_atoms = length(coords)
    for i in 2:(n_atoms - 1)
        angle_rad = Molly.bond_angle(coords[i-1], coords[i], coords[i+1], boundary)
        angle_deg = angle_rad * (180/pi)
        println("Bond angle at atom $i: $angle_deg degrees")
        push!(bond_angles,angle_deg)
    end
    return bond_angles
end

# Properties to differentiate the elements: 
# Mass, charge
n_atoms = 10
atom_mass = 15.9994u"g/mol"
boundary = CubicBoundary(1.5u"nm") 

atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=.2u"kJ * mol^-1") for i in 1:n_atoms] 
# sigma: distance where the intermolecular potential = 0. 
# or van der waals radius. how close two nonbonding particles can get. Minimum distance till molecules start messing with each other
# epsilon: well depth. how strongly they attract each other

temp = 100.0u"K"

coords = place_atoms(n_atoms ÷ 2, boundary; min_dist=0.3u"nm")
for i in 1:length(coords)
    push!(coords, coords[i] .+ [0.1, 0.0, 0.0]u"nm")
end

velocities = [random_velocity(atom_mass, temp) for i in 1:n_atoms]

#for bonds, if i want to calculate the bond angle i should use molly.bond_angle
bonds = InteractionList2Atoms(
    collect(1:(n_atoms ÷ 2)),           # First atom indices
    collect((1 + n_atoms ÷ 2):n_atoms), # Second atom indices
    [HarmonicBond(k=300_000.0u"kJ * mol^-1 * nm^-2", r0=0.1u"nm") for i in 1:(n_atoms ÷ 2)],
)

specific_inter_lists = (bonds,)

# All pairs apart from bonded pairs are eligible for non-bonded interactions
eligible = trues(n_atoms, n_atoms)
for i in 1:(n_atoms ÷ 2)
    eligible[i, i + (n_atoms ÷ 2)] = false
    eligible[i + (n_atoms ÷ 2), i] = false
end
neighbor_finder = DistanceNeighborFinder(
    eligible=eligible,
    n_steps=10,
    dist_cutoff=1.5u"nm",
)

cutoff = DistanceCutoff(1.2u"nm")
pairwise_inters = (LennardJones(use_neighbors=true, cutoff=cutoff),)

sys = System(
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    velocities=velocities,
    pairwise_inters=pairwise_inters,
    specific_inter_lists=specific_inter_lists,
    neighbor_finder=neighbor_finder,
    loggers=(
        temp=TemperatureLogger(10),
        coords=CoordinateLogger(5),

    ),
)

simulator = VelocityVerlet(
    dt=0.002u"ps",
    coupling=AndersenThermostat(temp, 1.0u"ps"),
)

simulate!(sys, simulator, 1_000)

neighbors = find_neighbors(sys)
show(neighbors)

visualize(
    sys.loggers.coords,
    boundary,
    "sim_diatomic.mp4";
    connections=[(i, i + (n_atoms ÷ 2)) for i in 1:(n_atoms ÷ 2)],
)