import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Timer;

public class Simulator {

    private final double deltaT;
    private double t=0.0;
    private final double maxT;
    private final ArrayList<Particle> particles;
    private final Scheme scheme;

    // Cantidad deseada de outputs totales (ajustable)
    private static final int targetOutputs = 150;

    private final double writeInterval;
    private double nextWriteTime;


    public Simulator(ArrayList<Particle> particles, Scheme scheme, double deltaT, double maxT, Path file) throws IOException {
        this.particles = particles;
        this.maxT = maxT;
        this.deltaT = deltaT;
        this.scheme = scheme;

        writeInterval = maxT / targetOutputs;
        nextWriteTime = writeInterval;
        executeSimulation(file);
    }

    public void executeSimulation(Path outputPath) throws IOException {
        try (OutputWriter out = OutputWriter.open(outputPath)) {

            // Para Beeman se necesita aceleración previa
            List<double[]> prevAcc = new ArrayList<>();
            // Para verlet original las posiciones iniciales
            List<double[]> prevPos = new ArrayList<>();

            switch (scheme) {
                case BEEMAN :
                    for (Particle p : particles) {
                        double[] f = Integrator.computeForce(p, particles);
                        prevAcc.add(new double[]{f[0]/p.getMass(), f[1]/p.getMass(), f[2]/p.getMass()});
                    };
                    break;

                case ORIGINAL_VERLET:
                    Integrator.initPrevPositions(particles, t, prevPos);
                    break;

                case GEAR_PREDICTOR_CORRECTOR_ORDER_5:
                    break;
            }

            while (t <= maxT) {
                if (t >= nextWriteTime) {
                    out.writeStep(particles, t);
                    nextWriteTime += writeInterval;
                }

                switch (scheme) {
                    case ORIGINAL_VERLET -> Integrator.updateParticlesVerlet(particles, deltaT, prevPos);
                    case BEEMAN -> Integrator.updateParticlesBeeman(particles, deltaT, prevAcc);
                    case GEAR_PREDICTOR_CORRECTOR_ORDER_5 -> Integrator.updateParticlesGear5(particles, deltaT);
                }

                t += deltaT;
                printProgress(t, maxT);
            }

        }
    }

    private void printProgress(double step, double maxT) {
        int percent = (int) ((step * 100.0) / maxT);
        String bar = "=".repeat(percent / 2) + " ".repeat(50 - percent / 2);
        System.out.printf("\r[%s] %d%%", bar, percent);
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 8) {
            throw new IllegalArgumentException("Wrong number of arguments");
        }
        int N = Integer.parseInt(args[0]);
        int iterations = Integer.parseInt(args[1]);
        String inputDir = args[2];
        String outputDir = args[3];
        Scheme scheme = Scheme.parseScheme(args[4]);
        String mode = args[5];
        double deltaT = Double.parseDouble(args[6]);
        int maxT = Integer.parseInt(args[7]);
        if (N <= 0 || iterations <= 0 || deltaT <=0 || maxT <= 0 || inputDir.isEmpty() || outputDir.isEmpty() || mode.isEmpty()) {
            System.out.print("Error: Parameters should be: N, iterations, inputDir, outputDir, Scheme(gp5/bm/vt), mode(single/cluster), deltaT, maxT\n");
            System.out.print("Example: 200 1 inputs outputs vt cluster 0.1 10");
            return;
        }
        Instant start = Instant.now();

//        while(N<=2000) {
//        for(double dt = deltaT; dt <= 0.09; dt *= 10) {
            for (int i = 0; i < iterations; i++) {
                InputParser parser = new InputParser(inputDir + "/" + mode + "/N" + N + "/input_N" + N + "_" + String.format("%04d", i) + ".csv", N);
                ArrayList<Particle> particles = parser.parseInputs();

                Path directory = Files.createDirectories(Path.of(outputDir, mode, Scheme.printScheme(scheme), "N_" + N));
                Path fileName = Path.of(directory + String.format("/output_N%d_%s_%s_t%d_%s.csv", N, mode, String.format(Locale.US, "dt%.3f", deltaT), maxT, String.format("%04d", i)));
                System.out.printf("\nStarting iteration %d/%d...\n", i + 1, iterations);
//                Simulator s = new Simulator(particles, scheme, deltaT, maxT, fileName);
                Simulator s = new Simulator(particles, scheme, deltaT, maxT, fileName);
                System.out.printf("\nIteration " + (i + 1) + " completed.");
            }
//        }
//            N+=100;
//        }

        Instant end = Instant.now();
        Duration elapsed = Duration.between(start, end);

        long hours = elapsed.toHours();
        long minutes = elapsed.toMinutesPart();
        long seconds = elapsed.toSecondsPart();

        System.out.printf("Tiempo transcurrido: %02d:%02d:%02d%n", hours, minutes, seconds);



    }



}
