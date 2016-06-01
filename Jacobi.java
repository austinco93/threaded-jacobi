import java.io.*;
import java.util.*;
/* Author: Austin Corotan
 * Date: 05/31/2016
 * Description: This program performs a threaded Jacobi Iteration using a barrier implemented with a monitor
 */
 
public class Jacobi {

    private static int matrixSize = 2048;
    private static double epsilon = 0.0001;
    
    /* This function reads in the input mtx file and returns a 2d array
     */
    public static double[][] readMTXFile(String mtxFileName) throws IOException {
    	BufferedReader stream = null;
    	int totalRows = matrixSize, totalColumns = matrixSize;
    	double[][] matrix = new double[totalRows][totalColumns];

    	try {
        	stream = new BufferedReader(new FileReader(mtxFileName));
        	for (int currentRow = 0; currentRow < totalRows; currentRow++) {
            	String line = stream.readLine();
            	String[] split = line.split(" ");
            	for (int currentColumn = 0; currentColumn < totalColumns; currentColumn++) {
                	matrix[currentRow][currentColumn] = Double.parseDouble(split[currentColumn]);
            	}
        	}
    	} finally {
        	if (stream != null) {
            	stream.close();
        	}
    	}
    	return matrix;
	}
        
    /* This function writes a 2d array to an mtx file 
     */
    public static void writeMTXFile(double[][] matrix) throws IOException {
    	BufferedWriter bw = null;
    	int totalRows = matrixSize, totalColumns = matrixSize;
    	
    	try {
        	File outputFile = new File("myoutput.mtx");
      		if (!outputFile.exists()) {
	     		outputFile.createNewFile();
	  		} 
	  		FileWriter fw = new FileWriter(outputFile);
        	bw = new BufferedWriter(fw);
        	for (int currentRow = 0; currentRow < totalRows; currentRow++) {
            	for (int currentColumn = 0; currentColumn < totalColumns; currentColumn++) {
                	bw.write(String.format("%.10f ",matrix[currentRow][currentColumn]));
            	}
            	bw.newLine();
        	}
    	} finally {
        	if (bw != null) {
            	bw.close();
        	}
    	}
    }
    
    /* this function will check if two 2d arrays are within a given threshold value (epsilon)
     */
    public static boolean verifyMTXFile(double[][] output, double[][] matrix, double epsilon) throws IOException {
        boolean result = true;
        for(int i = 1; i <= matrixSize-2; i++){
            for(int j = 1; j <= matrixSize-2; j++){
                if(Math.abs(output[i][j] - matrix[i][j]) < epsilon){
                    continue;
                }else{
                    System.out.printf("Failed at (%d,%d)\n", i, j);
                    result = false;
                    return result;
                }
            }
        }
        return result;
    }

    public static void main(String args[]) throws IOException, InterruptedException {
        int numThreads = 1;

        if(args.length == 1){
            numThreads = Integer.parseInt(args[0]);
        }

        ArrayList<Thread> threads = new ArrayList<Thread>(numThreads);
        double[][] jacMatrix = readMTXFile("input.mtx");
        double[][] jacMatrixCopy = readMTXFile("input.mtx");
        double[] maxDiff = new double[numThreads];
        int jacNumThreads = numThreads;
        int threadID = 1;
        double realMaxDiff = 0.001;

        Barrier barrier = new Barrier(numThreads);

        for(double max : maxDiff){
            max = 0.0;
        }
        long start = System.nanoTime();
        for (int x = 0; x < numThreads; x++) {
            Thread thread = new Thread(new JacobiThread(matrixSize, jacMatrix, jacMatrixCopy, jacNumThreads, maxDiff, epsilon, barrier, threadID, realMaxDiff));
            thread.start();
            threads.add(thread);
            threadID++;
        }
        for (Thread thread : threads) {
            thread.join();
        }
        long end = System.nanoTime();
        
        double[][] output = readMTXFile("output.mtx");

        if(verifyMTXFile(output, jacMatrix, epsilon)){
            System.out.printf("Completed Successfully!\n");
            System.out.printf("time = %d\n", end - start);
            writeMTXFile(jacMatrix);
        }		
    }
}

class JacobiThread implements Runnable {
    private double matrix[][];
    private double myCopy[][];
    private int numThreads;
    private double[] maxDiff;
    private double realMaxDiff;
    private double epsilon;
    private Barrier barrier = null;
    private int threadID;
    private int firstRow;
    private int lastRow;
    private int matrixSize;

    public JacobiThread(int matrixSize, double[][] matrix, double[][] myCopy, int numThreads, double[] maxDiff, double epsilon, Barrier barrier, int threadID, double realMaxDiff) {
        this.matrix = matrix;
        this.numThreads= numThreads;
        this.maxDiff = maxDiff;
        this.epsilon = epsilon;
        this.barrier = barrier;
        this.threadID = threadID;
        this.matrixSize = matrixSize;
        this.realMaxDiff = realMaxDiff;
        this.myCopy = myCopy;
    }
    
    @Override
    public void run(){
        try {
            firstRow = ((threadID - 1)*(matrixSize-2)/numThreads) + 1;
            lastRow = (threadID*(matrixSize-2))/numThreads;
            while(realMaxDiff > epsilon){
                for(int i = firstRow; i <= lastRow; i++){
                    for(int j = 1; j <= matrixSize-2; j++){
                        myCopy[i][j] = (matrix[i-1][j]+matrix[i+1][j]+matrix[i][j-1]+matrix[i][j+1])*0.25;
                    }
                }
                this.barrier.await();
                for(int i = firstRow; i <= lastRow; i++){
                    for(int j = 1; j <= matrixSize-2; j++){
                        matrix[i][j] = (myCopy[i-1][j]+myCopy[i+1][j]+myCopy[i][j-1]+myCopy[i][j+1])*0.25;
                    }
                }
                this.barrier.await();
                realMaxDiff = 0.0;
                maxDiff[threadID-1] = 0.0;
                for(int i = firstRow; i <= lastRow; i++){
                    for(int j = 1; j <= matrixSize-2; j++){
                        maxDiff[threadID-1] = Math.max(maxDiff[threadID-1], Math.abs(myCopy[i][j] - matrix[i][j]));
                    }
                }
                this.barrier.await();
                for(double max : maxDiff){
                    realMaxDiff = Math.max(realMaxDiff, max);
                }
            }
        }catch (InterruptedException ex){ 
            return; 
        }
    }
}

class Barrier {
    private int numThreads;
    private int awaitingThreads;

    public Barrier(int numThreads) {
        this.numThreads = numThreads;
        this.awaitingThreads = numThreads;
    }
    
    public synchronized void await() throws InterruptedException{
        awaitingThreads--;

        if(awaitingThreads > 0){
            this.wait();
        }else{
            awaitingThreads = numThreads;
            notifyAll();
        }
    }
}