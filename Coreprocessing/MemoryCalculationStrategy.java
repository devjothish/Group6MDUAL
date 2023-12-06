package Coreprocessing;

import java.util.concurrent.*;
import java.util.concurrent.atomic.DoubleAdder;
import java.util.logging.Level;
import java.util.logging.Logger;

// Memory Calculation Strategy
interface MemoryCalculationStrategy {
    double calculateMemoryUsage();
}

class DefaultMemoryCalculationStrategy implements MemoryCalculationStrategy {
    private static final double BYTES_PER_MB = 1024 * 1024;

    @Override
    public double calculateMemoryUsage() {
        Runtime runtime = Runtime.getRuntime();
        runtime.gc();
        return (double) (runtime.totalMemory() - runtime.freeMemory()) / BYTES_PER_MB;
    }
}

// Observer for Memory Events
interface MemoryObserver {
    void onPeakMemoryChanged(double newPeakMemory);
}

// Memory Monitor
public class MemoryMonitor {
    private final DoubleAdder maxMemory = new DoubleAdder();
    private final MemoryCalculationStrategy memoryStrategy;
    private final MemoryObserver observer;

    public MemoryMonitor(MemoryCalculationStrategy strategy, MemoryObserver observer) {
        this.memoryStrategy = strategy;
        this.observer = observer;
    }

    private void checkMemory() {
        double currentMemoryUsage = memoryStrategy.calculateMemoryUsage();
        double oldPeak = maxMemory.sum();
        if (currentMemoryUsage > oldPeak) {
            maxMemory.reset();
            maxMemory.add(currentMemoryUsage);
            observer.onPeakMemoryChanged(currentMemoryUsage);
        }
    }

    public void startMonitoring(long period, TimeUnit unit) {
        ScheduledExecutorService executor = Executors.newSingleThreadScheduledExecutor();
        executor.scheduleAtFixedRate(this::checkMemory, 0, period, unit);
    }

    public double getMaxMemoryUsage() {
        return maxMemory.sum();
    }
}
