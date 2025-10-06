import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GalaxyIntegrator {

    static private final double G = 1.0;
    static private final double H = 0.05;

    private static final double C0 = 3.0 / 20.0;
    private static final double C1 = 251.0 / 360.0;
    private static final double C2 = 1.0;
    private static final double C3 = 11.0 / 18.0;
    private static final double C4 = 1.0 / 6.0;
    private static final double C5 = 1.0 / 60.0;


    private GalaxyIntegrator() {
        throw new RuntimeException(this + ": Not instantiable");
    }

  static public List<double[]> computeForces(List<Particle> particles) {
      List<double[]> forces = new ArrayList<>();
      List<List<double[]>> interParticleForces = new ArrayList<>();
      for (int i = 0; i < particles.size() - 1; ++i) {
        Particle pi = particles.get(i);
        List<double[]> piForces = new ArrayList<>();

        for (int j = i + 1; j < particles.size(); ++j) {
          Particle pj = particles.get(j);

          double dx = pi.getX() - pj.getX();
          double dy = pi.getY() - pj.getY();
          double dz = pi.getZ() - pj.getZ();

          double dist2 = dx * dx + dy * dy + dz * dz + H * H;
          double dist3 = Math.pow(dist2, 1.5);

          double factor = -G * pi.getMass() * pj.getMass() / dist3;

          piForces.add(new double[] { factor * dx, factor * dy, factor * dz });
        }
        interParticleForces.add(piForces);
      }
      // [
      // p0: [p1, p2, p3, p4, p5],
      // p1: [p2, p3, p4, p5],
      // p2: [p3, p4, p5],
      // p3: [p4, p5],
      // p4: [p5],
      // ]
      for (int i = 0; i < particles.size(); ++i) {
        var totalForce = new double[] { 0, 0, 0 };
        if (i < particles.size() - 1) {
          var forces1 = interParticleForces.get(i);
          for (var force : forces1) {
            totalForce[0] += force[0];
            totalForce[1] += force[1];
            totalForce[2] += force[2];
          }
        }
        for (int j = i - 1; j >= 0; --j) {
          var forces2 = interParticleForces.get(j);
          var force = forces2.get(i - j - 1);
          totalForce[0] -= force[0];
          totalForce[1] -= force[1];
          totalForce[2] -= force[2];
        }
        forces.add(totalForce);
      }

      return forces;
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

    static public void updateParticlesBeeman(List<Particle> particles, double dt,
        // Paso de integración con Beeman
        List<double[]> prevAccelerations) {
      // Aceleraciones actuales
      List<double[]> accNow = new ArrayList<>();
      var forces = computeForces(particles);
      for (int i = 0; i < particles.size(); ++i) {
        var f = forces.get(i);
        var p = particles.get(i);
        accNow.add(new double[] { f[0] / p.getMass(), f[1] / p.getMass(), f[2] / p.getMass() });
      }

      // Actualizar posiciones
      for (int i = 0; i < particles.size(); i++) {
        Particle p = particles.get(i);
        double[] aNow = accNow.get(i);
        double[] aPrev = prevAccelerations.get(i);

        double newX = p.getX() + p.getVx() * dt + (2.0 / 3.0) * aNow[0] * dt * dt - (1.0 / 6.0) * aPrev[0] * dt * dt;
        double newY = p.getY() + p.getVy() * dt + (2.0 / 3.0) * aNow[1] * dt * dt - (1.0 / 6.0) * aPrev[1] * dt * dt;
        double newZ = p.getZ() + p.getVz() * dt + (2.0 / 3.0) * aNow[2] * dt * dt - (1.0 / 6.0) * aPrev[2] * dt * dt;

        p.setPosition(newX, newY, newZ);
      }

      // Aceleraciones futuras (con posiciones nuevas)
      List<double[]> accNext = new ArrayList<>();
      forces = computeForces(particles);
      for (int i = 0; i < particles.size(); ++i) {
        var f = forces.get(i);
        var p = particles.get(i);
        accNext.add(new double[] { f[0] / p.getMass(), f[1] / p.getMass(), f[2] / p.getMass() });
      }

      // Actualizar velocidades
      for (int i = 0; i < particles.size(); i++) {
        Particle p = particles.get(i);
        double[] aPrev = prevAccelerations.get(i);
        double[] aNow = accNow.get(i);
        double[] aNext = accNext.get(i);

        double newVx = p.getVx() + (1.0 / 3.0) * aNext[0] * dt + (5.0 / 6.0) * aNow[0] * dt - (1.0 / 6.0) * aPrev[0] * dt;
        double newVy = p.getVy() + (1.0 / 3.0) * aNext[1] * dt + (5.0 / 6.0) * aNow[1] * dt - (1.0 / 6.0) * aPrev[1] * dt;
        double newVz = p.getVz() + (1.0 / 3.0) * aNext[2] * dt + (5.0 / 6.0) * aNow[2] * dt - (1.0 / 6.0) * aPrev[2] * dt;

        p.setVelocity(newVx, newVy, newVz);
      }

      // Importante: devolver accNow como "prev" para el próximo paso
      prevAccelerations.clear();
      prevAccelerations.addAll(accNow);
      //
    }


    // Velocity-Verlet
    public static void updateParticlesVelocityVerlet(List<Particle> particles, double dt) {

        // 1. Calcular fuerzas actuales
        List<double[]> fNow = computeForces(particles);

        // 2. Actualizar posiciones en el lugar
        for (int i = 0; i < particles.size(); i++) {
            Particle p = particles.get(i);
            double[] f = fNow.get(i);
            double ax = f[0] / p.getMass();
            double ay = f[1] / p.getMass();
            double az = f[2] / p.getMass();

            p.setPosition(
                    p.getX() + p.getVx() * dt + 0.5 * ax * dt * dt,
                    p.getY() + p.getVy() * dt + 0.5 * ay * dt * dt,
                    p.getZ() + p.getVz() * dt + 0.5 * az * dt * dt
            );
        }

        // 3. Calcular nuevas fuerzas con posiciones actualizadas
        List<double[]> fNext = computeForces(particles);

        // 4. Actualizar velocidades
        for (int i = 0; i < particles.size(); i++) {
            Particle p = particles.get(i);
            double[] f0 = fNow.get(i);
            double[] f1 = fNext.get(i);

            p.setVelocity(
                    p.getVx() + 0.5 * dt * (f0[0] + f1[0]) / p.getMass(),
                    p.getVy() + 0.5 * dt * (f0[1] + f1[1]) / p.getMass(),
                    p.getVz() + 0.5 * dt * (f0[2] + f1[2]) / p.getMass()
            );
        }
    }


    public static void initGear(List<Particle> particles, double dt, double[][][] r, double[][][] p) {
        int N = particles.size();
        r[0] = new double[6][N]; r[1] = new double[6][N]; r[2] = new double[6][N];
        p[0] = new double[6][N]; p[1] = new double[6][N]; p[2] = new double[6][N];
        for (int i = 0; i < N; i++) {
            Particle part = particles.get(i);
            r[0][0][i] = part.getX();  r[1][0][i] = part.getY();  r[2][0][i] = part.getZ();
            r[0][1][i] = part.getVx(); r[1][1][i] = part.getVy(); r[2][1][i] = part.getVz();
        }
        for (int i = 0; i < N; i++) {
            Particle pi = particles.get(i);
            double[] f = computeForce(pi, particles);
            r[0][2][i] = f[0] / pi.getMass();
            r[1][2][i] = f[1] / pi.getMass();
            r[2][2][i] = f[2] / pi.getMass();
        }
        for (int k = 3; k <= 5; k++) {
            Arrays.fill(r[0][k], 0.0);
            Arrays.fill(r[1][k], 0.0);
            Arrays.fill(r[2][k], 0.0);
        }
    }

    public static void updateParticlesGear5(List<Particle> particles, double dt, double[][][] r, double[][][] p) {
        final int N = particles.size();
        final double dt1 = dt, dt2 = dt1*dt, dt3 = dt2*dt, dt4 = dt3*dt, dt5 = dt4*dt;
        for (int i = 0; i < N; i++) {
            // X
            double r0 = r[0][0][i], r1 = r[0][1][i], r2 = r[0][2][i], r3 = r[0][3][i], r4 = r[0][4][i], r5 = r[0][5][i];
            p[0][0][i] = r0 + r1*dt1 + r2*dt2/2.0 + r3*dt3/6.0 + r4*dt4/24.0 + r5*dt5/120.0;
            p[0][1][i] = r1 + r2*dt1 + r3*dt2/2.0 + r4*dt3/6.0 + r5*dt4/24.0;
            p[0][2][i] = r2 + r3*dt1 + r4*dt2/2.0 + r5*dt3/6.0;
            p[0][3][i] = r3 + r4*dt1 + r5*dt2/2.0;
            p[0][4][i] = r4 + r5*dt1;
            p[0][5][i] = r5;

            // Y
            r0=r[1][0][i]; r1=r[1][1][i]; r2=r[1][2][i]; r3=r[1][3][i]; r4=r[1][4][i]; r5=r[1][5][i];
            p[1][0][i] = r0 + r1*dt1 + r2*dt2/2.0 + r3*dt3/6.0 + r4*dt4/24.0 + r5*dt5/120.0;
            p[1][1][i] = r1 + r2*dt1 + r3*dt2/2.0 + r4*dt3/6.0 + r5*dt4/24.0;
            p[1][2][i] = r2 + r3*dt1 + r4*dt2/2.0 + r5*dt3/6.0;
            p[1][3][i] = r3 + r4*dt1 + r5*dt2/2.0;
            p[1][4][i] = r4 + r5*dt1;
            p[1][5][i] = r5;

            // Z
            r0=r[2][0][i]; r1=r[2][1][i]; r2=r[2][2][i]; r3=r[2][3][i]; r4=r[2][4][i]; r5=r[2][5][i];
            p[2][0][i] = r0 + r1*dt1 + r2*dt2/2.0 + r3*dt3/6.0 + r4*dt4/24.0 + r5*dt5/120.0;
            p[2][1][i] = r1 + r2*dt1 + r3*dt2/2.0 + r4*dt3/6.0 + r5*dt4/24.0;
            p[2][2][i] = r2 + r3*dt1 + r4*dt2/2.0 + r5*dt3/6.0;
            p[2][3][i] = r3 + r4*dt1 + r5*dt2/2.0;
            p[2][4][i] = r4 + r5*dt1;
            p[2][5][i] = r5;
        }

        for (int i = 0; i < N; i++) {
            particles.get(i).setPosition(p[0][0][i], p[1][0][i], p[2][0][i]);
        }

        double[] R2x = new double[N], R2y = new double[N], R2z = new double[N];
        for (int i = 0; i < N; i++) {
            Particle pi = particles.get(i);
            double[] fP = computeForce(pi, particles);
            double axP = fP[0] / pi.getMass(), ayP = fP[1] / pi.getMass(), azP = fP[2]/ pi.getMass();;

            double dAx = axP - p[0][2][i];
            double dAy = ayP - p[1][2][i];
            double dAz = azP - p[2][2][i];

            R2x[i] = dAx * dt2 * 0.5;
            R2y[i] = dAy * dt2 * 0.5;
            R2z[i] = dAz * dt2 * 0.5;
        }

        final double invDt  = 1.0 / dt;
        final double invDt2 = invDt * invDt;
        final double invDt3 = invDt2 * invDt;
        final double invDt4 = invDt3 * invDt;
        final double invDt5 = invDt4 * invDt;

        for (int i = 0; i < N; i++) {
            // X
            double r0C = p[0][0][i] + C0 * R2x[i];
            double r1C = p[0][1][i] + C1 * (R2x[i] * invDt);
            double r2C = p[0][2][i] + C2 * (2.0 * R2x[i] * invDt2);
            double r3C = p[0][3][i] + C3 * (6.0 * R2x[i] * invDt3);
            double r4C = p[0][4][i] + C4 * (24.0 * R2x[i] * invDt4);
            double r5C = p[0][5][i] + C5 * (120.0 * R2x[i] * invDt5);

            r[0][0][i]=r0C; r[0][1][i]=r1C; r[0][2][i]=r2C; r[0][3][i]=r3C; r[0][4][i]=r4C; r[0][5][i]=r5C;

            // Y
            r0C = p[1][0][i] + C0 * R2y[i];
            r1C = p[1][1][i] + C1 * (R2y[i] * invDt);
            r2C = p[1][2][i] + C2 * (2.0 * R2y[i] * invDt2);
            r3C = p[1][3][i] + C3 * (6.0 * R2y[i] * invDt3);
            r4C = p[1][4][i] + C4 * (24.0 * R2y[i] * invDt3 * invDt);
            r5C = p[1][5][i] + C5 * (120.0 * R2y[i] * invDt3 * invDt2);

            r[1][0][i]=r0C; r[1][1][i]=r1C; r[1][2][i]=r2C; r[1][3][i]=r3C; r[1][4][i]=r4C; r[1][5][i]=r5C;

            // Z
            r0C = p[2][0][i] + C0 * R2z[i];
            r1C = p[2][1][i] + C1 * (R2z[i] * invDt);
            r2C = p[2][2][i] + C2 * (2.0 * R2z[i] * invDt2);
            r3C = p[2][3][i] + C3 * (6.0 * R2z[i] * invDt3);
            r4C = p[2][4][i] + C4 * (24.0 * R2z[i] * invDt3 * invDt);
            r5C = p[2][5][i] + C5 * (120.0 * R2z[i] * invDt3 * invDt2);

            r[2][0][i]=r0C; r[2][1][i]=r1C; r[2][2][i]=r2C; r[2][3][i]=r3C; r[2][4][i]=r4C; r[2][5][i]=r5C;
        }

        for (int i = 0; i < N; i++) {
            particles.get(i).setPosition(r[0][0][i], r[1][0][i], r[2][0][i]);
            particles.get(i).setVelocity(r[0][1][i], r[1][1][i], r[2][1][i]);
        }
    }



//    public static void initPrevPositions(List<Particle> particles, double dt, List<double []> prev) {
//
//        for (Particle p : particles) {
//            double[] f = computeForce(p, particles);
//            double ax = f[0] / p.getMass();
//            double ay = f[1] / p.getMass();
//            double az = f[2] / p.getMass();
//
//            double xPrev = p.getX() - p.getVx()*dt + 0.5*ax*dt*dt;
//            double yPrev = p.getY() - p.getVy()*dt + 0.5*ay*dt*dt;
//            double zPrev = p.getZ() - p.getVz()*dt + 0.5*az*dt*dt;
//
//            prev.add(new double[]{xPrev, yPrev, zPrev});
//        }
//    }


    // Un paso de Velocity-Verlet para la partícula pi
//    static public void updateParticle(Particle pi, List<Particle> particles, double dt) {
//
//        double[] fOld = computeForce(pi, particles);
//
//        double newX = pi.getX() + pi.getVx() * dt + 0.5 * dt * dt * fOld[0] / pi.getMass();
//        double newY = pi.getY() + pi.getVy() * dt + 0.5 * dt * dt * fOld[1] / pi.getMass();
//        double newZ = pi.getZ() + pi.getVz() * dt + 0.5 * dt * dt * fOld[2] / pi.getMass();
//
//        // Crear "fantasma" para calcular nuevas fuerzas
//        Particle temp = new Particle(pi.getId(), newX, newY, newZ, pi.getVx(), pi.getVy(), pi.getVz(), pi.getMass());
//        double[] fNew = computeForce(temp, particles);
//
//        double newVx = pi.getVx() + 0.5 * (fOld[0] + fNew[0]) / pi.getMass() * dt;
//        double newVy = pi.getVy() + 0.5 * (fOld[1] + fNew[1]) / pi.getMass() * dt;
//        double newVz = pi.getVz() + 0.5 * (fOld[2] + fNew[2]) / pi.getMass() * dt;
//
//        pi.setPosition(newX, newY, newZ);
//        pi.setVelocity(newVx, newVy, newVz);
//    }


//    static public void updateParticles(List<Particle> particles, double dt) {
//
//        for (Particle p : particles) {
//            updateParticle(p, particles, dt);
////            updateParticleFree(p, dt);
//        }
//    }

    // Verlet original
//    public static void updateParticlesVerlet(List<Particle> particles, double dt,
//                                             List<double[]> prevPositions) {
//
//        List<double[]> fNow = new ArrayList<>();
//
//        for (Particle p : particles) {
//            fNow.add(computeForce(p, particles));
//        }
//
//        List<double[]> posNow = new ArrayList<>();
//
//        for (int i = 0; i < particles.size(); i++) {
//
//            Particle p = particles.get(i);
//            posNow.add(new double[]{p.getX(), p.getY(), p.getZ()});
//
//            double[] f = fNow.get(i);
//
//            double[] prevPosition = prevPositions.get(i);
//
//            double newX = 2*p.getX() - prevPosition[0] + (dt*dt / p.getMass())*f[0];
//            double newY = 2*p.getY() - prevPosition[1] + (dt*dt / p.getMass())*f[1];
//            double newZ = 2*p.getZ() - prevPosition[2] + (dt*dt / p.getMass())*f[2];
//
//            double newVx = (newX-prevPosition[0])/(2*dt);
//            double newVy = (newY-prevPosition[1])/(2*dt);
//            double newVz = (newZ-prevPosition[2])/(2*dt);
//
//            p.setPosition(newX, newY, newZ);
//            p.setVelocity(newVx, newVy, newVz);
//        }
//
//        prevPositions.clear();
//        prevPositions.addAll(posNow);
//    }

}

