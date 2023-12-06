package Coreprocessing;

import java.lang.management.*;

public interface PerformanceMetric {
    long measure();
}

class CPUPerformanceMetric implements PerformanceMetric {
    @Override
    public long measure() {
        ThreadMXBean bean = ThreadMXBeanFactory.getInstance();
        return bean.isCurrentThreadCpuTimeSupported() ? bean.getCurrentThreadCpuTime() : 0L;
    }
}

class MemoryPerformanceMetric implements PerformanceMetric {
    private static final int BYTES_PER_MB = 1024 * 1024;
    private static MemoryPerformanceMetric instance = new MemoryPerformanceMetric();
    private long peakMemoryUsage = 0;

    private MemoryPerformanceMetric() {}

    public static MemoryPerformanceMetric getInstance() {
        return instance;
    }

    @Override
    public long measure() {
        long currentUsage = calculateMemoryUsage();
        updatePeakMemoryUsage(currentUsage);
        return currentUsage / BYTES_PER_MB;
    }

    private long calculateMemoryUsage() {
        Runtime runtime = Runtime.getRuntime();
        runtime.gc();
        return runtime.totalMemory() - runtime.freeMemory();
    }

    private void updatePeakMemoryUsage(long currentUsage) {
        if (currentUsage > peakMemoryUsage) {
            peakMemoryUsage = currentUsage;
        }
    }

    public long getPeakMemoryUsage() {
        return peakMemoryUsage / BYTES_PER_MB;
    }
}

class ThreadMXBeanFactory {
    private static ThreadMXBean beanInstance = ManagementFactory.getThreadMXBean();

    private ThreadMXBeanFactory() {}

    public static ThreadMXBean getInstance() {
        return beanInstance;
    }
}
