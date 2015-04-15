
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.Semaphore;

public class Jacobi {

	private int N;
	private int numProcs;
	private double LEFT = 1.0;
	private double TOP = 1.0;
	private double RIGHT = 80.0;
	private double BOTTOM = 80.0;
	private double EPSILON = 0.10;
	private double[][] grid;
	private double[][] newGrid;
	private int HEIGHT;
	private int iters = 0;
	private int MAXITERS = 20000;
	private boolean debug = true;
	private double maxDiff[];
	private int numArrived = 0;
	private Semaphore barrier, e;
	private long startTime;
	private long endTime;
	private long endMicroSeconds;
	private long endSeconds;
	private FileWriter out = null;
	
	public static void main(String[] args) throws IOException {
		new Jacobi(args);
	}
	
	public Jacobi(String[] args) throws IOException
	{
		if(args.length == 0)
		{
			usage();
			return;
		}
		else
		{
			setN(Integer.parseInt(args[0]));
			setNumProcs(Integer.parseInt(args[1]));
		}
		if(args.length > 2)
		{
			if(args.length != 7)
			{
				usage();
				return;
			}
			if(args[2] != null)
				setLEFT(Double.parseDouble(args[2]));
			if(args[3] != null)
				setTOP(Double.parseDouble(args[3]));
			if(args[4] != null)
				setRIGHT(Double.parseDouble(args[4]));
			if(args[5] != null)
				setBOTTOM(Double.parseDouble(args[5]));
			if(args[6] != null)
				setEPSILON(Double.parseDouble(args[6]));
		}
		setMaxDiff();
		setHeight();
		setGrid();
		
		
		barrier = new Semaphore(1);
		e = new Semaphore(numProcs);
		
		JacobiWorker[] threads = new JacobiWorker[numProcs];
		
		startTime = System.nanoTime();
		for(int i = 0; i < numProcs; i ++)
		{
			threads[i] = new JacobiWorker(i);
			threads[i].start();
		}
//		endTime = System.nanoTime();
		
		for(int i = 0; i < numProcs; i ++)
		{
			try {
				threads[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			
		}
		endTime = System.nanoTime();
		endSeconds = (endTime - startTime)/1000000000;
		endMicroSeconds = ((endTime - startTime) - (endSeconds * 1000000000))/1000;
		
		try {
			out = new FileWriter("JacobiResults.txt");
			out.write("Grid\t = " + N + "x" + N + "\n");
			out.write("numProcs = " + numProcs + "\n");
			out.write("left\t = " + LEFT + "\n");
			out.write("top\t\t = " + TOP + "\n");
			out.write("right\t = " + RIGHT + "\n");
			out.write("bottom\t = " + BOTTOM + "\n");
			out.write("epsilon\t = " + EPSILON + "\n");
			out.write("execution time = " + endSeconds + " seconds, " + endMicroSeconds + " microseconds" + "\n\n");
			out.write(retGrid(grid));
		} finally {
			if(out != null) out.close();
		}
		
		System.out.println("main: numProcs = " + numProcs + ", N = " + N);
		System.out.println("execution time = " + endSeconds + " seconds, " + endMicroSeconds + " microseconds" + "\n\n");
	}
	
	private void barrier(int id)
	{
		//lock
		try {
			e.acquire();
			barrier.acquire();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		numArrived++;
		if(numArrived == numProcs)
		{
			numArrived = 0;
			//start
			e.release(numProcs);
		} else
			try {
				Thread.sleep(50);
			} catch (InterruptedException e1) {
				e1.printStackTrace();
			}//wait
		barrier.release();
		//unlock
	}
	
	
	public static void usage()
	{
		System.out.println("JacobiJava Usage\n");
		System.out.println("JacobiJava N numProcs [L] [T] [R] [B] [E]\n");
		System.out.println("  N is the size of the grid that will be NxN (required)");
		System.out.println("  numProcs is the number of processes (threads) to create to solve the problem (required, bounded by minimum of 1, maximum of 16)");
		System.out.println("  L is the fixed value for the left edge of the grid [optional, default = 1.0]");
		System.out.println("  T is the fixed value for the top edge of the grid [optional, default = 1.0]");
		System.out.println("  R is the fixed value for the right edge of the grid [optional, default = 80.0]");
		System.out.println("  B is the fixed value for the bottom edge of the grid [optional, default = 80.0]");
		System.out.println("  E is the epsilon value used to determine when to stop the computation [optional, default = 1.0]\n");
		System.out.println("This program will print out the settings to a file called JacobiJava.log");
		System.out.println("as well as the final grid contents.\n");
		System.out.println("The number of processes, iterations, and computation time will be printed out");
		System.out.println("to standard output.\n");
	}
	
	public void setGrid()
	{
		grid = new double[N+2][N+2];
		newGrid = new double[N+2][N+2];
		for(int i = 0; i < N + 2; i ++)
		{
			for(int j = 0; j < N + 2; j ++)
			{
				if(i == 0 && j != 0)
					grid[i][j] = TOP;
				else if(j == 0)
					grid[i][j] = LEFT;
				else if(i == N + 1)
					grid[i][j] = BOTTOM;
				else if(j == N + 1)
					grid[i][j] = RIGHT;
				else
					grid[i][j] = 0.0;
				newGrid[i][j] = grid[i][j];
			}
		}
		
		if(!debug)
		{
			printGrid(grid);
			System.out.println();
			printGrid(newGrid);
		}
	}
	
	public void printMaxDiff()
	{
		for(int i = 0; i < numProcs; i ++)
		{
			System.out.print(maxDiff[i] + " ");
		}
		System.out.println();
	}
	
	public void printGrid(double[][] grid)
	{
		for(int i = 0; i < N + 2; i ++)
		{
			for(int j = 0; j < N + 2; j ++)
			{
				System.out.printf("%2.4f ", grid[i][j]);
			}
			System.out.println();
		}
	}
	
	public String retGrid(double[][] grid)
	{
		String gridString = "";
		for(int i = 0; i < N + 2; i ++)
		{
			gridString += "\t";
			for(int j = 0; j < N + 2; j ++)
			{
				
				gridString += String.format("%05.6g", grid[i][j]) + " ";
			}
			gridString += "\n";
		}
		return gridString;
	}

	public int getN() {
		return N;
	}

	public int getNumProcs() {
		return numProcs;
	}

	public double getLEFT() {
		return LEFT;
	}

	public double getTOP() {
		return TOP;
	}

	public double getRIGHT() {
		return RIGHT;
	}

	public double getBOTTOM() {
		return BOTTOM;
	}

	public double getEPSILON() {
		return EPSILON;
	}
	
	public int getHeight()
	{
		return HEIGHT;
	}

	public void setN(int n) {
		N = n;
	}

	public void setNumProcs(int numProcs) {
		this.numProcs = numProcs;
	}

	public void setLEFT(double lEFT) {
		LEFT = lEFT;
	}

	public void setTOP(double tOP) {
		TOP = tOP;
	}

	public void setRIGHT(double rIGHT) {
		RIGHT = rIGHT;
	}

	public void setBOTTOM(double bOTTOM) {
		BOTTOM = bOTTOM;
	}

	public void setEPSILON(double ePSILON) {
		EPSILON = ePSILON;
	}
	
	public void setHeight()
	{
		HEIGHT = N/numProcs;
	}
	
	public void setMaxDiff()
	{
		maxDiff = new double[numProcs];
		for(int i = 0; i < numProcs; i ++)
		{
			maxDiff[i] = 0.0;
		}
	}

	class JacobiWorker extends Thread {
			
			private int proc;
			
			public JacobiWorker(int proc)
			{
				this.proc = proc;
			}
			
			public void run()
			{
				int firstRow = (proc) * HEIGHT + 1;
				int lastRow = firstRow + HEIGHT ;
				double mydiff = 0.0;
				
				barrier(proc);
				for(iters = 1; iters < MAXITERS; iters++)
				{
					for(int i = firstRow; i < lastRow; i ++)
					{
						for(int j = 1; j < N + 1; j ++)
						{
							newGrid[i][j] = (grid[i-1][j] + grid[i+1][j] + grid[i][j-1] + grid[i][j+1]) * .25;
						}
					}
					
					barrier(proc);
					
					for(int i = firstRow; i < lastRow; i ++)
					{
						for(int j = 1; j < N + 1; j ++)
						{
							grid[i][j] = (newGrid[i-1][j] + newGrid[i+1][j] + newGrid[i][j-1] + newGrid[i][j+1]) * .25;
							
							mydiff = Math.max(mydiff, Math.abs(grid[i][j] - newGrid[i][j]));
						}
					}
					maxDiff[proc] = mydiff;
					if(mydiff < EPSILON)
						break;
					mydiff = 0.0;
					barrier(proc);
				}
				
				if(!debug)
				{
					printGrid(grid);
					System.out.println();
					printMaxDiff();
				}
			}
		}
}
