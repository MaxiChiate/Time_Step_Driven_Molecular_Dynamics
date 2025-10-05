public enum Scheme {
    GEAR_PREDICTOR_CORRECTOR_ORDER_5,
    BEEMAN,
    VERLET;

    public static Scheme parseScheme(String s) {
        if (s == null) throw new IllegalArgumentException("scheme is null");
        return switch (s.trim().toLowerCase()) {
            case "gp5" -> GEAR_PREDICTOR_CORRECTOR_ORDER_5;
            case "bm"  -> BEEMAN;
            case "vt"  -> VERLET;
            default -> throw new IllegalArgumentException("Unknown scheme: " + s);
        };
    }

    public static String printScheme(Scheme scheme){
        return switch (scheme){
            case GEAR_PREDICTOR_CORRECTOR_ORDER_5 ->  "gear_order_5";
            case BEEMAN ->  "beeman";
            case VERLET ->   "verlet";
        };
    }

}
