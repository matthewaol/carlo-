using Molly
using GLMakie
using Statistics
# Properties to differentiate the elements: 
n_atoms = 100
atom_mass = 15.9994u"g/mol" # oxygen mass
boundary = CubicBoundary(2.0u"nm") 
atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=.2u"kJ * mol^-1") for i in 1:n_atoms] 

# sigma: distance where the intermolecular potential = 0. 
# or van der waals radius. how close two nonbonding particles can get. Minimum distance till molecules start messing with each other
# epsilon: well depth. how strongly they attract each other

temp = 10000u"K"

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
        coords=CoordinateLogger(10),
        press =PressureLogger(10)
    ),
)

simulator = VelocityVerlet(
    dt=0.002u"ps",
    coupling=AndersenThermostat(temp, 1.0u"ps"),
    #BerendsenThermostat(temp, .002u"ps"),
    #AndersenThermostat(temp, 1.0u"ps")

)

simulate!(sys, simulator, 6_0000)

temps = ustrip(values(sys.loggers.temp))
pressures = ustrip(values(sys.loggers.press))
mean_temp = mean(temps)
mean_pressures = mean(pressures)

pv = pressures * (2.0*10^-9)^2
nrt = fill(50*.08314, length(pv)) .* temps 

println("Mean temp is: ", mean_temp)
println("Mean pressure is: ", mean_pressures)

f = Figure(resolution=(600, 400))
ax = Axis(
    f[1, 1],
    xlabel="nrt",
    ylabel="pv",
    title="pv vs. nrt",
)
scatter!(
    nrt,
    pv,
    markersize=5,
)
save("pv_nrt.png", f)

#temp fig
f = Figure(resolution=(600, 400))
ax = Axis(
    f[1, 1],
    xlabel="Step",
    ylabel="Temperature",
    title="Temperature change during sim",
)
scatter!(
    ax,
    100 .* (1:length(temps)),
    temps,
    markersize=5,
)
save("temp_oxygen.png", f)

#press fig
f = Figure(resolution=(600, 400))
ax = Axis(
    f[1, 1],
    xlabel="Step",
    ylabel="Pressure",
    title="Pressure change during sim",
)
scatter!(
    ax,
    100 .* (1:length(pressures)),
    pressures,
    markersize=5,
)
save("press_oxygen.png", f)


visualize(
    sys.loggers.coords,
    boundary,
    "sim_diatomic.mp4";
    connections=[(i, i + (n_atoms ÷ 2)) for i in 1:(n_atoms ÷ 2)],
)