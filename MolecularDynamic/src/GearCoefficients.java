public class GearCoefficients {

    public static final double C0_GALAXY = 3.0 / 20.0;
    public static final double C0_OSCILLATOR = 3.0/16.0;
    public static final double C1 = 251.0 / 360.0;
    public static final double C2 = 1.0;
    public static final double C3 = 11.0 / 18.0;
    public static final double C4 = 1.0 / 6.0;
    public static final double C5 = 1.0 / 60.0;

    private GearCoefficients() {
        throw new RuntimeException(this.getClass().getSimpleName() + " cannot be instantiated.");
    }
}
