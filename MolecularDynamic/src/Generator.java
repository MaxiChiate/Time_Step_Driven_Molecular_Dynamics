import java.io.IOException;
import java.util.Locale;
import java.util.Random;

public class Generator {
    private static final double MASS = 1.0;    // todas las partículas
    // Defaults
    private static final double velocity = 0.1; // v0 o vx0
    //For clusters building
    private static final double dx = 4.0;
    private static final double dy = 0.5;

    private int particleId = 0;

    private final int N;
    private final Mode mode;
    private final int iterations;
    private final double v0;     // single
    private final double vx0;    // clusters

    public enum Mode { SINGLE, CLUSTERS }

    public Generator(int N, Mode mode, int iterations, double v0_or_vx0) {
        this.N = N;
        this.mode = mode;
        this.iterations = iterations;
        if (mode == Mode.SINGLE) {
            this.v0 = v0_or_vx0;
            this.vx0 = Double.NaN;
        } else {
            this.v0 = Double.NaN;
            this.vx0 = v0_or_vx0;
        }

    }

    //Elegir la dirección de la velocidad es tomar un punto aleatorio en la superficie de la espera imaginaria
    //de radio |v|. No se pueden hacer 3 randoms directamente para armar el vector porque los puntos se van a distribuir más
    //fuerte cerca de los polos cuando hago el random z
    //https://stackoverflow.com/questions/14805583/dispersing-n-points-uniformly-on-a-sphere
    //https://mathworld.wolfram.com/SpherePointPicking.html
    static double[] randomVelocity(Random rng, double speed) {
        double u = 2.0 * rng.nextDouble() - 1.0;       // u uniforme en [-1,1]
        double phi = 2.0 * Math.PI * rng.nextDouble(); // phi uniforme en [0, 2pi]
        double s = Math.sqrt(1.0 - u*u);
        return new double[] {
                speed * s * Math.cos(phi),
                speed * s * Math.sin(phi),
                speed * u };
    }

    //La idea es que si tomo 3 valores INDEPENDIENTES que tengan Dist Normal con media=0 y desvío=1 =>
    //la posición cumple con tener media 0 y desvio estandar 1
    private void generateSingle(Random rng, int iter) throws IOException {
        Particle[] particles = new Particle[N];
        particleId = 0;
        for (int i = 0; i < N; i++) {
            double x = rng.nextGaussian();
            double y = rng.nextGaussian();
            double z = rng.nextGaussian();
            double[] u = randomVelocity(rng, v0);
            double vx = u[0];
            double vy = u[1];
            double vz = u[2];
            particles[i] = new Particle(particleId++,x, y, z, vx, vy, vz, MASS);
        }

        OutputWriter.writeOutputForSingle(iter, particles, N);
    }

    private void generateClusters(Random rng, int iter) throws IOException {
        Particle[] particleClusterA = new Particle[N];
        Particle[] particleClusterB = new Particle[N];
        particleId = 0;
        // Cada cúmulo ~ N(0,1) alrededor de su centro:
        // 1: (-dx/2, -dy/2, 0), 2: (+dx/2, +dy/2, 0)
        for (int i = 0; i < N; i++) {
            double x = rng.nextGaussian() - dx / 2.0;
            double y = rng.nextGaussian() - dy / 2.0;
            double z = rng.nextGaussian();
            particleClusterA[i] = new Particle(particleId++,x, y, z, vx0, 0, 0, MASS);
        }
        for (int i = 0; i < N; i++) {
            double x = rng.nextGaussian() + dx / 2.0;
            double y = rng.nextGaussian() + dy / 2.0;
            double z = rng.nextGaussian();
            particleClusterB[i] = new Particle(particleId++,x, y, z, -vx0, 0, 0, MASS);
        }
        OutputWriter.writeOutputForClusters(iter, particleClusterA, particleClusterB, N);
    }

    public void generateAll() throws IOException {
        for (int i = 0; i < iterations; i++) {
            Random rng = new Random();
            if (mode == Mode.SINGLE) {
                generateSingle(rng, i);
            } else {
                generateClusters(rng, i);
            }
        }
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 3) {
            System.out.print("Usage: N iterations single/clusters\n");
            return;
        }
        int N = Integer.parseInt(args[0]);
        int iterations = Integer.parseInt(args[1]);
        Mode mode = switch (args[2].toLowerCase(Locale.ROOT)) {
            case "single" -> Mode.SINGLE;
            case "clusters" -> Mode.CLUSTERS;
            default -> throw new IllegalArgumentException();
        };
        Generator gen = new Generator(N, mode, iterations, velocity);
        gen.generateAll();
    }
}
