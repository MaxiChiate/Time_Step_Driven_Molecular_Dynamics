import java.util.ArrayList;
import java.util.List;

public class Integrator {

    static private final double G = 1.0;
    static private final double H = 0.05;

    private Integrator() {
        throw new RuntimeException(this + ": Not instantiable");
    }

    static public double[] computeForce(Particle pi, List<Particle> particles) {
        double fx = 0, fy = 0, fz = 0;

        for (Particle pj : particles) {
            if (pi.getId() == pj.getId()) continue;

            double dx = pi.getX() - pj.getX();
            double dy = pi.getY() - pj.getY();
            double dz = pi.getZ() - pj.getZ();

            double dist2 = dx * dx + dy * dy + dz * dz + H * H;
            double dist3 = Math.pow(dist2, 1.5);

            double factor = -G * pi.getMass() * pj.getMass() / dist3;

            fx += factor * dx;
            fy += factor * dy;
            fz += factor * dz;
        }

        return new double[]{fx, fy, fz};
    }

    // Un paso de Velocity-Verlet para la partícula pi
    static public void updateParticle(Particle pi, List<Particle> particles, double dt) {

        double[] fOld = computeForce(pi, particles);

        double newX = pi.getX() + pi.getVx() * dt + 0.5 * dt * dt * fOld[0] / pi.getMass();
        double newY = pi.getY() + pi.getVy() * dt + 0.5 * dt * dt * fOld[1] / pi.getMass();
        double newZ = pi.getZ() + pi.getVz() * dt + 0.5 * dt * dt * fOld[2] / pi.getMass();

        // Crear "fantasma" para calcular nuevas fuerzas
        Particle temp = new Particle(pi.getId(), newX, newY, newZ, pi.getVx(), pi.getVy(), pi.getVz(), pi.getMass());
        double[] fNew = computeForce(temp, particles);

        double newVx = pi.getVx() + 0.5 * (fOld[0] + fNew[0]) / pi.getMass() * dt;
        double newVy = pi.getVy() + 0.5 * (fOld[1] + fNew[1]) / pi.getMass() * dt;
        double newVz = pi.getVz() + 0.5 * (fOld[2] + fNew[2]) / pi.getMass() * dt;

        pi.setPosition(newX, newY, newZ);
        pi.setVelocity(newVx, newVy, newVz);
    }

    // Paso de integración con Beeman
    static public void updateParticlesBeeman(List<Particle> particles, double dt,
                                             List<double[]> prevAccelerations) {
        // Aceleraciones actuales
        List<double[]> accNow = new ArrayList<>();
        for (Particle p : particles) {
            double[] f = computeForce(p, particles);
            accNow.add(new double[]{f[0]/p.getMass(), f[1]/p.getMass(), f[2]/p.getMass()});
        }

        // Actualizar posiciones
        for (int i = 0; i < particles.size(); i++) {
            Particle p = particles.get(i);
            double[] aNow = accNow.get(i);
            double[] aPrev = prevAccelerations.get(i);

            double newX = p.getX() + p.getVx()*dt + (2.0/3.0)*aNow[0]*dt*dt - (1.0/6.0)*aPrev[0]*dt*dt;
            double newY = p.getY() + p.getVy()*dt + (2.0/3.0)*aNow[1]*dt*dt - (1.0/6.0)*aPrev[1]*dt*dt;
            double newZ = p.getZ() + p.getVz()*dt + (2.0/3.0)*aNow[2]*dt*dt - (1.0/6.0)*aPrev[2]*dt*dt;

            p.setPosition(newX, newY, newZ);
        }

        // Aceleraciones futuras (con posiciones nuevas)
        List<double[]> accNext = new ArrayList<>();
        for (Particle p : particles) {
            double[] f = computeForce(p, particles);
            accNext.add(new double[]{f[0]/p.getMass(), f[1]/p.getMass(), f[2]/p.getMass()});
        }

        // Actualizar velocidades
        for (int i = 0; i < particles.size(); i++) {
            Particle p = particles.get(i);
            double[] aPrev = prevAccelerations.get(i);
            double[] aNow = accNow.get(i);
            double[] aNext = accNext.get(i);

            double newVx = p.getVx() + (1.0/3.0)*aNext[0]*dt + (5.0/6.0)*aNow[0]*dt - (1.0/6.0)*aPrev[0]*dt;
            double newVy = p.getVy() + (1.0/3.0)*aNext[1]*dt + (5.0/6.0)*aNow[1]*dt - (1.0/6.0)*aPrev[1]*dt;
            double newVz = p.getVz() + (1.0/3.0)*aNext[2]*dt + (5.0/6.0)*aNow[2]*dt - (1.0/6.0)*aPrev[2]*dt;

            p.setVelocity(newVx, newVy, newVz);
        }

        // Importante: devolver accNow como "prev" para el próximo paso
        prevAccelerations.clear();
        prevAccelerations.addAll(accNow);
    }


    static public void updateParticles(List<Particle> particles, double dt) {

        for (Particle p : particles) {
            updateParticle(p, particles, dt);
//            updateParticleFree(p, dt);
        }
    }
}