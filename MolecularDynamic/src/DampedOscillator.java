import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.Format;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.util.Locale;

public class DampedOscillator {
    private final static double k = 10000;
    private final static double m = 70;
    private final static double gamma = 100;
    private final double dt;
    private final static double maxT = 5;

    private double r0, r1, r2, r3, r4, r5;
    private double r0P, r1P, r2P, r3P, r4P, r5P;

    private double x;
    private double v;
    private double t;
    private double aPrev;
    private double xPrev;
    private double vPrev;
    private final double printDt = 0.01;
    private double tToPrint = 0.0;

    public double acceleration (double x, double v){
        return (-k*x - gamma*v)/m;
    }

    private void initBeeman(double x0, double v0){
        this.aPrev = acceleration(x0, v0);
    }
    //Con Euler Modificado
    private void initVerlet(double x0, double v0){
        this.vPrev = v0 - dt * acceleration(x0, v0);
        this.xPrev = x0 - v0 * dt + 0.5 * acceleration(x,v0) * dt * dt;
    }

    public DampedOscillator(double x0, double v0, double dt) {
        this.dt=dt;
        this.x=x0;
        this.v=v0;
        this.t=0;
    }

    private void stepBeeman() {
        double aNow = acceleration(x, v);

        double xNext = x + v * dt + (2.0/3.0) * aNow * dt*dt - (1.0/6.0) * aPrev * dt*dt;

        double vPred = v + 1.5 * aNow * dt - 0.5 * aPrev * dt;

        double aNext = acceleration(xNext, vPred);

        double vNext = v + (1.0/3.0) * aNext * dt + (5.0/6.0) * aNow * dt - (1.0/6.0) * aPrev * dt;

        x = xNext;
        v = vNext;
        aPrev = aNow;
        t += dt;
    }

    private void stepOriginalVerlet() throws IOException {
        //Como necesito la fuerza y esta depende de la velocidad uso velocidad anterior como estimativo (v siempre vale VPrev excpeto para t=0)
        double aNow  = acceleration(x, v);
        double xNext = 2.0 * x - xPrev + aNow * dt * dt;
        double vNow = (xNext - xPrev) / (2.0 * dt);

        vPrev = vNow;
        xPrev = x;
        x     = xNext;
        v     = vNow;
        t    += dt;
    }



    private void predictOrder5(){
        double dt1 = dt;
        double dt2 = dt1*dt;
        double dt3 = dt2*dt;
        double dt4 = dt3*dt;
        double dt5 = dt4*dt;

        r0P = r0 + r1*dt1 + r2*dt2/2.0 + r3*dt3/6.0 + r4*dt4/24.0 + r5*dt5/120.0;
        r1P = r1 + r2*dt1 + r3*dt2/2.0 + r4*dt3/6.0 + r5*dt4/24.0;
        r2P = r2 + r3*dt1 + r4*dt2/2.0 + r5*dt3/6.0;
        r3P = r3 + r4*dt1 + r5*dt2/2.0;
        r4P = r4 + r5*dt1;
        r5P = r5;
    }

    private void initGear(double x0, double v0) {
        r0=x0;
        r1=v0;
        r2=acceleration(x0,v0);
        r3 = -(k/m) * r1 - (gamma/m) * r2;
        r4 = -(k/m) * r2 - (gamma/m) * r3;
        r5 = -(k/m) * r3 - (gamma/m) * r4;
    }
    private void stepGearOrder5() {

        final double C0 = GearCoefficients.C0_OSCILLATOR;
        final double C1 = GearCoefficients.C1;
        final double C2 = GearCoefficients.C2;
        final double C3 = GearCoefficients.C3;
        final double C4 = GearCoefficients.C4;
        final double C5 = GearCoefficients.C5;

        predictOrder5();
        double aP = acceleration(r0P, r1P);
        double a = aP - r2P;

        double R2 = a*dt*dt/2;

        double r0C=r0P + C0 * R2;
        double r1C=r1P + C1 * R2/dt;
        double r2C=r2P + C2 * R2*2.0/Math.pow(dt,2);
        double r3C=r3P + C3 * R2*6.0/Math.pow(dt,3);
        double r4C=r4P + C4 * R2*24.0/Math.pow(dt,4);
        double r5C=r5P + C5 * R2*120.0/Math.pow(dt,5);

        x = r0C;
        v = r1C;
        t += dt;

        r0=r0C;
        r1=r1C;
        r2=r2C;
        r3=r3C;
        r4=r4C;
        r5=r5C;
    }

    private void stepTheoretical(){
        t+=dt;
        double beta = gamma/(2*m);
        double wd = Math.pow((k / m) - beta*beta, 0.5);
        x = 1 * Math.exp(-beta*t)*Math.cos(wd*t);
        v = 1 * Math.exp(-beta*t)*(-beta * Math.cos(wd *t) - wd * Math.sin(wd *t));
    }

    public double x()        { return x; }
    public double vx()       { return v; }
    public double t()        { return t; }
    public double deltaT()   { return dt; }
    public double maxT()     { return maxT; }
    public double mass()     { return m; }
    public double vPrev()    { return vPrev;}
    public double xPrev()    { return xPrev;}

    public void simulate(Path directory, Scheme scheme) throws IOException {
        Files.createDirectories(directory);
        String dtString = String.format(Locale.US,"%.9f", dt);
        String fileName = switch (scheme) {

            case BEEMAN -> "beeman"+dtString+".csv";
            case GEAR_PREDICTOR_CORRECTOR_ORDER_5 -> "gear_order_5"+dtString+".csv";
            case VERLET -> "original_verlet"+dtString+".csv";
            case null -> "theoretical"+dtString+".csv";
        };
        try (OutputWriter out = OutputWriter.open(directory.resolve(fileName))) {
            while (t() <= maxT()) {
                switch (scheme) {
                    case BEEMAN -> {
                        if (t() >= tToPrint) {
                            out.writeStepDamped(this, t());
                            tToPrint += printDt;
                        }
                        stepBeeman();
                    }
                    case GEAR_PREDICTOR_CORRECTOR_ORDER_5 -> {
                        if (t() >= tToPrint) {
                            out.writeStepDamped(this, t());
                            tToPrint += printDt;
                        }
                        stepGearOrder5();
                    }
                    case VERLET -> {
                        if (Double.compare(t(), 0) >= 0) {
                            stepOriginalVerlet();
                            if (t() - dt >= tToPrint) {
                                out.writeStepDampedVerlet(this, t() - dt);
                                tToPrint += printDt;
                            }
                        } else {
                            stepOriginalVerlet();
                        }
                    }
                    case null -> {
                        if (t() >= tToPrint) {
                            out.writeStepDamped(this, t());
                            tToPrint += printDt;
                        }
                        stepTheoretical();
                    }
                }
            }
        }
    }

    public static void main(String[] args) throws IOException {
        double x0 = 1.0;
        double v0 = - gamma/(2*m);

        //simulates until 10^-8
        int maxPower = 8;

        for (int e = 1; e <= maxPower; e++) {
            double dt1 = Math.pow(10.0, -e);
            if (e != 1) {
                simulateForDT(x0, v0, dt1);
            }
            if (e < maxPower) {
                double dt2 = 5.0 * Math.pow(10.0, -(e + 1));
                simulateForDT(x0, v0, dt2);
            }
        }
    }

    private static void simulateForDT(double x0, double v0, double dt) throws IOException {
        DampedOscillator oscBeeman = new DampedOscillator(x0, v0, dt);
        oscBeeman.initBeeman(x0, v0);
        DampedOscillator oscGear = new DampedOscillator(x0, v0, dt);
        oscGear.initGear(x0, v0);
        DampedOscillator oscVerlet = new DampedOscillator(x0, v0, dt);
        oscVerlet.initVerlet(x0, v0);
        DampedOscillator oscTheo = new DampedOscillator(x0, v0, dt);

        Path dir = Path.of("outputs/dampedOscillator");
        oscBeeman.simulate(dir, Scheme.BEEMAN);
        oscGear.simulate(dir, Scheme.GEAR_PREDICTOR_CORRECTOR_ORDER_5);
        oscVerlet.simulate(dir, Scheme.VERLET);
        oscTheo.simulate(dir, null);
    }
}