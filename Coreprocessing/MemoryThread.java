package Coreprocessing;
import java.util.logging.Level;
import java.util.logging.Logger;

public class MemoryThread extends Thread {
	private static final double MegaBytes = 1024*1024;
    public double maxMemory = 0;
    public void computeMemory() {
        Runtime.getRuntime().gc();
        double used = (Runtime.getRuntime().totalMemory()- Runtime.getRuntime().freeMemory())/MegaBytes;
        if(maxMemory < used)
            maxMemory = used;
    }

    @Override
    public void run() {
        while (true) {
            computeMemory();
            try {
                Thread.sleep(100);
            } catch (InterruptedException ex) {
                Logger.getLogger(MemoryThread.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

}
