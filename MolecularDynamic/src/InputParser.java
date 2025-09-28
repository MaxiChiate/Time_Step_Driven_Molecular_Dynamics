import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import static java.lang.System.exit;

public class InputParser {
    private final String inputPath;
    private final int particles_amount;

    public InputParser(String inputPath, int particles_amount) {
        this.inputPath = inputPath;
        this.particles_amount = particles_amount;
    }
    public ArrayList<Particle> parseInputs() {
        ArrayList<Particle> particles = new ArrayList<>();
        Path path = Path.of(inputPath);
        if(!Files.exists(path)) {
            throw new  IllegalArgumentException("File not found: " + inputPath);
        }
        try {
            List<String> lines = Files.readAllLines(path);
            for (String line : lines) {
                String[] parts = line.split(",");
                if (parts.length != 8) {
                    throw new IllegalArgumentException("Expected: columns " +
                            "(id,x,y,z,vx,vy,vz,mass) and found: " + parts.length +
                            ". Content: '" + line + "'");
                }
                try {
                    Particle particle = getParticle(parts);
                    particles.add(particle);
                } catch (NumberFormatException e) {
                    System.out.printf("Number Format exception: '%s'\n", line);
                }
            }
            if (particles.size() != particles_amount) {
                throw new IllegalArgumentException("Number of particles does not match the expected amount");
            }
        } catch (Exception e) {
            System.out.print("Error while parsing input\n");
            exit(1);
        }
        return  particles;
    }

    private static Particle getParticle(String[] parts) {
        int id = Integer.parseInt(parts[0]);
        double x = Double.parseDouble(parts[1]);
        double y = Double.parseDouble(parts[2]);
        double z = Double.parseDouble(parts[3]);
        double vx = Double.parseDouble(parts[4]);
        double vy = Double.parseDouble(parts[5]);
        double vz = Double.parseDouble(parts[6]);
        double mass = Double.parseDouble(parts[7]);
        Particle particle = new Particle(id, x, y, z, vx, vy, vz, mass);
        return particle;
    }
}