public class Particle {
    private int id;
    private double x;
    private double y;
    private double z;
    private double vx;
    private double vy;
    private double vz;
    private double mass;

    public Particle(int id, double x, double y, double z, double vx, double vy, double vz, double mass) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.z = z;
        this.vx = vx;
        this.vy = vy;
        this.vz = vz;
        this.mass = mass;
    }

    public void setPosition(double x, double y, double z) {
        this.x = x; this.y = y; this.z = z;
    }

    public void setVelocity(double vx, double vy, double vz) {
        this.vx = vx; this.vy = vy; this.vz = vz;
    }

    public int getId() {
        return id;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getZ() {
        return z;
    }

    public double getVx() {
        return vx;
    }

    public double getVy() {
        return vy;
    }

    public double getVz() {
        return vz;
    }

    public double getMass() {
        return mass;
    }
}
