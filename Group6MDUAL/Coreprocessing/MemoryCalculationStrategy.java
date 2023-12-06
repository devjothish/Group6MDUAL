package Coreprocessing;

import java.util.*;
import java.util.logging.*;

public class MemoryCalculationStrategy extends Thread {
    private static final double BYTES_PER_MB = 1024 * 1024;
    public static int megaMB = 0;
    private double maxMemoryUsage = 0;

    private void checkAndUpdatePeakMemoryUsage() {
        double currentUsage = getCurrentMemoryUsage();
        maxMemoryUsage = Math.max(maxMemoryUsage, currentUsage);
    }

    private double getCurrentMemoryUsage() {
        Runtime runtime = Runtime.getRuntime();
        runtime.gc();
        return (runtime.totalMemory() - runtime.freeMemory()) / BYTES_PER_MB;
    }

    @Override
    public void run() {
        while (!isInterrupted()) {
            checkAndUpdatePeakMemoryUsage();
            pauseMonitoring();
        }
    }

    private void pauseMonitoring() {
        try {
            Thread.sleep(100);
        } catch (InterruptedException ex) {
            Logger.getLogger(MemoryCalculationStrategy.class.getName()).log(Level.SEVERE, "Thread interrupted", ex);
            Thread.currentThread().interrupt();
        }

    }


    public double getPeakMemoryUsage() {
        return maxMemoryUsage;
    }

    public void beginMonitoring() {
    }
}