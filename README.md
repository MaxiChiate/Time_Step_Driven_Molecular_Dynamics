# Gravitational Molecular Dynamics

## Objective
Simulate the gravitational collapse and interaction of N-body systems (galaxies or clusters) using time-driven molecular dynamics and analyze the stability and energy conservation of different integration schemes.

## Methodology
- Time-driven simulation (fixed timestep $\Delta t$)
- Gravitational force calculation with smoothing parameter $H$ to avoid singularities: $F = -G \frac{m_i m_j}{(r^2 + H^2)^{3/2}} \mathbf{r}$
- Integration schemes: **Gear Predictor-Corrector (Order 5)**, **Beeman**, and **Velocity Verlet**
- Parallelized force computation using Java Streams and ForkJoinPool
- Two initial configurations: **Single** (Gaussian distribution) and **Clusters** (two colliding systems)

## Results
- System transitions from initial collapse to a stationary regime defined by the half-mass radius $r_{hm}$
- Total energy is conserved within expected limits depending on the integration scheme and $\Delta t$
- Gear Predictor-Corrector (GP5) shows superior stability and lower energy error for small timesteps
- Time to reach equilibrium ($t^*$) and final distribution depend on particle number $N$ and initial velocity

## Conclusion
The simulator provides a robust framework for N-body gravitational simulations, validating the implementation of high-order integrators. Macroscopic evolution, such as the virialization of the system, emerges from microscopic gravitational interactions, with Gear Order 5 being the most precise scheme for long-term stability studies.

---

## Analysis

- **Presentation**  
  Full methodology, equations, plots, and quantitative analysis are documented in the slides used for the oral defense:  
  `analysis/SdS_TP4_2025Q2G02_Presentación.pdf`

- **Simulation Videos**
  - **Two-cluster collision — Gear Predictor-Corrector (Order 5)**  
    https://youtu.be/3yjdBzG9-yQ  

  The video shows the time evolution of two initially separated gravitational clusters, highlighting collapse, interaction, and post-collision relaxation under a high-order integration scheme.


---

## Execution (Reproducibility)

### Compile
```bash
javac MolecularDynamic/src/*.java -d out/production/SDS4
```

### Generate Initial Conditions
```bash
# Usage: java -cp out/production/SDS4 Generator <N> <iterations> <single/clusters>

# Example for single system:
java -cp out/production/SDS4 Generator 1000 1 single

# Example for two clusters:
java -cp out/production/SDS4 Generator 500 1 clusters
```

### Run Simulation
```bash
# Usage: java -cp out/production/SDS4 Simulator <N> <iterations> <inputDir> <outputDir> <scheme> <mode> <deltaT> <maxT> <targetOutputs>
# Schemes: gp5 (Gear), bm (Beeman), vt (Verlet)

# Example for single system:
java -cp out/production/SDS4 Simulator 1000 1 ./MolecularDynamic/inputs ./MolecularDynamic/outputs gp5 single 0.0001 50 500

# Example for two clusters:
java -cp out/production/SDS4 Simulator 500 1 ./MolecularDynamic/inputs ./MolecularDynamic/outputs gp5 clusters 0.0001 50 500
```

### Optional Visualization
```bash
# Requires ffmpeg for mp4 export

# Visualizing single system:
python3 visualizer.py MolecularDynamic/outputs/single/gear_order_5/N1000/output_N1000_single_dt0.00010_t50_0000.csv --fps 30 --export

# Visualizing clusters (uses --mode clusters for pink/blue coloring):
python3 visualizer.py MolecularDynamic/outputs/clusters/gear_order_5/N0500/output_N500_clusters_dt0.00010_t50_0000.csv --mode clusters --fps 30 --export
```
