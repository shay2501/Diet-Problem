import java.io.*;

/**
 *
 * @author Arijit Ghosh
 */
public class DietProblem {
    
    static int M;                   // no of constraints
    static int N;                   // no of variables
    static int c;
    static int d;
    static int[] basicVariable;     // basic variable corresponding to row i
    
    static double[][] a;
    static double[][] aTranspose;
    static double[] b;
    static double[] cj;
    /******************************************************************************************/
    public static void main(String args[]){
        
        long startTime = System.currentTimeMillis();
        
        String file = null;
         if(args.length > 0)
         {
             file = args[0];
         }
         
        double[] y = readInputs(file);
        
        if(y != null)
        {
            for (int i = 0; i < y.length; i++) 
            {
                System.out.println("y" + (i+1) + " = " + y[i]);
             }
        }
        //System.out.println("time to execute simplex = " + (System.currentTimeMillis() - startTime));
    }
    /******************************************************************************************/
    private static double[] readInputs(String file) {
            
        //return value for our solution
        double[] y = null;
        BufferedReader br = null;
        
        try{
            String line = null;
            // Reading the input file
            if(file != null)
            {
                br = new BufferedReader(new FileReader(file)); 
            }
            else
            {
                br = new BufferedReader(new FileReader("file100_100.txt"));                         
            }
            
            line = br.readLine();
            while(line.startsWith("//"))
            {
              line = br.readLine();
            }
            String[] tokens = line.split(" ");
            cj = new double[tokens.length];
            for(int i = 0; i < tokens.length; i++) {
                cj[i] = Double.parseDouble(tokens[i]);
               
            }
            N = cj.length;
           
            
            line = br.readLine();
            while(line.startsWith("//"))
            {
              line = br.readLine();
            }
            tokens = line.split(" ");    
            b = new double[tokens.length];
            for(int i = 0; i < tokens.length; i++) {
                b[i] = Double.parseDouble(tokens[i]);
             
            }
            M = b.length;
     
            

            a = new double[M+1][N+M+1];
            for (int i = 0; i < M; i++) {
                line = br.readLine();
                while(line.startsWith("//"))
                {
                  line = br.readLine();
                }
                tokens = line.split(" ");
                for (int j = 0; j < N; j++) {
                    a[i][j] = Double.parseDouble(tokens[j]);     
                    
                }     
            }
            
            aTranspose = new double[N+M+1][M+1];
            for ( c = 0 ; c < M ; c++ )
            {
               for ( d = 0 ; d < N ; d++ ) {
                    aTranspose[d][c] = a[c][d];                    
                }              
            }   
            
            for ( c = 0 ; c < N ; c++ )
            {
               for ( d = 0 ; d < M ; d++ ) {                    
                  
                }
            }           

            //create vectors for Basic and Non-Basic vars
            int[] NB = new int[b.length];
            int[] B = new int[cj.length];
            
            //call simplex and get solution to maximization of transpose
            double[] x = Simplex.simplex(NB, B, aTranspose, cj, b);
            
            //get the objective function from the result to calculate our minimum solution
            double[] objFn = Simplex.returnObjFn();
            
            //number of non basic vars in the primal objective function
            int n = objFn.length;
            
            //new array to hold the solution to the diet problem, the size of the non-basic array
            y = new double[N];
            
            //foreach solution variable, see if it is in the NB objective function
            for(int i = 0; i< N; i++)
            {
                //for every initial non basic variable, search the current non basic variable array
                for(int j = 0; j< NB.length; j++)
                {
                   //if the variable to calculate is in the non basic vector, its coef is negation of the 
                   //its coef in the objective function from the dual solution
                     if(NB[j] == (n+i+1))
                     {
                         y[i] = (-1 * objFn[j]);
                     }
                }
            }

        }
        catch(FileNotFoundException e){
            System.err.println("Input file missing or of different format: " + e.getMessage());
            
        }
        catch(IOException e){
            System.err.println("There was an error reading from the file: " + e.getMessage());
        }

        return y;
    }
    /******************************************************************************************/
}
    
   