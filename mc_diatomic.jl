# Properties to differentiate the elements: 
using Molly
using GLMakie
using Statistics
# Properties to differentiate the elements: 
n_atoms = 150
atom_mass = 15.9994u"g/mol" # oxygen mass
boundary = CubicBoundary(2.0u"nm") 
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
    [HarmonicBond(k=500_000.0u"kJ * mol^-1 * nm^-2", r0=0.1u"nm") for i in 1:(n_atoms ÷ 2)],
) # 500kJ for bond energy of oxygen

specific_inter_lists = (bonds,)

temperatures = [1198.0, 798.0, 398.0, 198.0, 98.0, 8.0]u"K"
sys = System(   
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    #pairwise_inters=pairwise_inters,
    specific_inter_lists = specific_inter_lists,
    loggers=(
        coords=CoordinateLogger(n_atoms, dims=n_dimensions(boundary)),
        montecarlo=MonteCarloLogger(),
        temp=TemperatureLogger(10),
        press=PressureLogger(10)
    ),
)

trial_args = Dict(:shift_size => 0.1u"nm")
for t in temperatures
    sim = MetropolisMonteCarlo(; 
        temperature=t,
        trial_moves=random_uniform_translation!,
        trial_args=trial_args,
    )

    simulate!(sys, sim, 10_00)
end

println(sys.loggers.montecarlo.n_accept)


visualize(
    sys.loggers.coords,
    boundary,
    "sim_montecarlo.mp4";
    connections=[(i, i + (n_atoms ÷ 2)) for i in 1:(n_atoms ÷ 2)],
)
