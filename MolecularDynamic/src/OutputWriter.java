import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Formatter;
import java.util.List;
import java.util.Locale;

public class OutputWriter implements AutoCloseable {

    private static final String INPUTS_PATH = "./inputs";
    private final BufferedWriter bw;
    private final StringBuilder sb;
    private final Formatter fmt;
//    private double tToPrint=0;
//    private static final double printDt=0.001;

    private OutputWriter(Path path) throws IOException {
        this.bw = Files.newBufferedWriter(
                path,
                StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING,
                StandardOpenOption.WRITE
        );
        this.sb  = new StringBuilder();
        this.fmt = new Formatter(sb, Locale.US);
    }

    public static OutputWriter open(Path path) throws IOException {
        return new OutputWriter(path);
    }

//--------FOR DAMPED OSCILATOR-------------------------------------------------------------

 public void writeStepDamped(DampedOscillator d, double time) throws IOException {
         sb.setLength(0);
         fmt.format("%.4f%n", time);
         fmt.format(
                 "%.17g,%.17g,%.5f%n",
                 d.x(), d.vx(), d.mass());
         bw.write(sb.toString());
 }
    public void writeStepDampedVerlet(DampedOscillator d, double t) throws IOException {
            sb.setLength(0);
            fmt.format("%.4f%n", t);
            fmt.format(
                    "%.17g,%.17g,%.5f%n",
                    d.xPrev(), d.vPrev(), d.mass());
            bw.write(sb.toString());
    }

//---------FOR GALAXY SIMULATION ----------------------------------------------------------
public void writeStep(List<Particle> particles, double time) throws IOException {
    sb.setLength(0);
    fmt.format("%.4f%n", time);
    for (Particle p : particles) {
        fmt.format("%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.5f%n", p.getId(), p.getX(), p.getY(), p.getZ(), p.getVx(), p.getVy(), p.getVz(), p.getMass());
    }
    bw.write(sb.toString());
}

public void writeAnalysis(double time, double energy, double rhm, Double tStar) throws IOException {
    sb.setLength(0);

    if (tStar == null)
        fmt.format(Locale.US, "%.6f,%.17g,%.17g,%s%n", time, energy, rhm, "NaN");
    else
        fmt.format(Locale.US, "%.6f,%.17g,%.17g,%.6f%n", time, energy, rhm, tStar);
    bw.write(sb.toString());
}



//----------FOR GALAXY INPUT GENERATION----------------------------------------------------
    public static void writeOutputForSingle(int iter, Particle[] particles, int N) throws IOException {
        File file = prepareOutputFile("single", iter, N);
        StringBuilder sb = new StringBuilder();
        appendParticlesFormatted(sb, particles, N);
        writeToFile(file, sb);
    }

    public static void writeOutputForClusters(int iter, Particle[] particles1, Particle[] particles2, int N) throws IOException {
        File file = prepareOutputFile("clusters", iter, N);
        StringBuilder sb = new StringBuilder();
        appendParticlesFormatted(sb, particles1, N); // mismo bucle con N
        appendParticlesFormatted(sb, particles2, N); // mismo bucle con N
        writeToFile(file, sb);
    }

    private static File prepareOutputFile(String modeFolder, int iter, int N) throws IOException {
        String dirPath = INPUTS_PATH + "/" + modeFolder + "/N" + String.format("%04d", N);
        Files.createDirectories(Path.of(dirPath));
        String fileName = String.format("input_N%d_%s.csv", N, String.format("%04d", iter));
        File file = new File(dirPath + "/" + fileName);
        if (file.exists()) file.delete();
        file.createNewFile();
        return file;
    }

    private static void writeToFile(File file, StringBuilder sb) throws IOException {
        Files.writeString(file.toPath(), sb.toString(), StandardOpenOption.TRUNCATE_EXISTING);
    }

    private static void appendParticlesFormatted(StringBuilder sb, Particle[] particles, int N) {
        for (int i = 0; i < N; i++) {
            sb.append(String.format(Locale.US,
                    "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.5g%n",
                    particles[i].getId(),
                    particles[i].getX(), particles[i].getY(), particles[i].getZ(),
                    particles[i].getVx(), particles[i].getVy(), particles[i].getVz(), particles[i].getMass()));
        }
    }

    @Override
    public void close() throws IOException {
        fmt.close();
        bw.close();
    }
}
