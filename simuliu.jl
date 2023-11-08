n_atoms = 100
atoms = [Atom(mass=10.0u"g/mol", charge=1.0) for i in 1:n_atoms]
boundary = RectangularBoundary(4.0u"nm")

coords = place_atoms(n_atoms, boundary; min_dist=0.2u"nm")
pairwise_inters = (Coulomb(),)

temperatures = [1198.0, 798.0, 398.0, 198.0, 98.0, 8.0]u"K"
sys = System(
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    pairwise_inters=pairwise_inters,
    loggers=(
        coords=CoordinateLogger(n_atoms, dims=n_dimensions(boundary)),
        montecarlo=MonteCarloLogger(),
    ),
)

trial_args = Dict(:shift_size => 0.1u"nm")
for t in temperatures
    sim = MetropolisMonteCarlo(; 
        temperature=t,
        trial_moves=random_uniform_translation!,
        trial_args=trial_args,
    )

    simulate!(sys, sim, 10_000)
end

println(sys.loggers.montecarlo.n_accept)
# 15234

visualize(sys.loggers.coords, boundary, "sim_montecarlo.gif")