
import java.io.BufferedReader;
import java.io.IOException;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.LinkedList;

/*
 * This source file implements the Simplex Algorithm for solving Linear Programming problems. 
 * It has a main function, and can be run on it's own with an input file. 
 * The simplex method expects and empty array for Basic and NonBasic variables, and a populated
 * primal A matrix consisting of coefficients of the constraints of the LP, a populated c vector
 * that holds the objective function coefficients and a populated b vector with the RHS of the constraints. 
 * 
 * It returns a solution vector for the initial non-basic variables
 */

/**
 *
 * @author Shannon Whalen
 */
public class Simplex {
    
    private static double[] objFn;
    
    /******************************************************************************************/
     /*
     * Main allows for the capability of taking an input file and reading it line by line
     * to test the algorithm independent of any application
     */
    public static void main(String args[]){
        //adjacency matrix to hold the graph
         double[][] A = null; //matrix
         double [] b = null; // vector for the LHS
         double [] c = null; //vector for the objective function
         
         long startTime = System.currentTimeMillis();

         String line = null;
       
         int lineNum = 0;
         
         String file = null;
         if(args.length > 0)
         {
             file = args[0];
         }

         try
         {
           BufferedReader f = null;
           if(file == null)
           {
            f = new BufferedReader(new InputStreamReader(System.in));
           }
           else
           {
             f = new BufferedReader(new FileReader(file));
           }
           
           //read in line by line
           line = f.readLine();
           
           int idx = 0;
           while(line != null && line.trim().length() > 0)
           {
               
               //if this is a comment, then read another line
               while(line != null && line.startsWith("//"))
               {
                   line = f.readLine();
               }
               
               String[] tokens = line.split(" ");
               
               if(lineNum == 0)
               {
                   //initialize c
                   c = new double[tokens.length];
               
                   //reading in the objective coefficient
                   for(int i = 0; i< tokens.length; i++)
                   {
                       c[i] = Double.parseDouble(tokens[i]);
                   }
               
               }
                   
               else if(lineNum == 1)
               {
                   //intialize b
                   b = new double[tokens.length];
                   
                   //reading in the objective coefficient
                   for(int i = 0; i< tokens.length; i++)
                   {
                       b[i] = Double.parseDouble(tokens[i]);
                   }
               }
               else{
              
                if(A == null && c != null && b != null)
                {
                    A = new double[b.length][c.length];
                }

                //populate the matrix with the input values
                  
                for(int i = 0; i< tokens.length; i++)
                {
                    A[idx][i] = Double.parseDouble(tokens[i]);
                }

                //read the next line and increment the row index
                idx++;
              }
               
              //read the next line
              line = f.readLine();
              
              //increase the line number
              lineNum++;
           }
           
           //return the basic and non-basic vars from init simplex
           int N[] = new int[c.length];
           int B[] = new int[b.length];
           
           //simplex returns the solution vector for the non-basic variables
           double[] x = simplex(N, B, A, b, c);
           
           //if not unbounded or infeasible, print out solution to non-basic vars
           if(x != null)
           {
               for(int i = 0; i< x.length; i++)
               {
                 System.out.println("x" + (i+1) + " = " + x[i]);
               }
           }
           
           //System.out.println("Time to execute Simplex: " + (System.currentTimeMillis() - startTime));
           
         }
         catch (IOException e) {
           System.out.println("Error: " + e.getMessage());
           System.exit(1);
         }
         
 }
    /******************************************************************************************/
    /*
     * returnObjFn returns the objective function vector in case a solution to the dual needs calculating
     */
    public static double[] returnObjFn()
    {
        //return the objective function
        return objFn;
    }
    /******************************************************************************************/
    /*
     * Simplex solves a linear programming problem defined by a coefficient matrix A, a constraint
     * vector b and an objective function c. Vectors N and B hold the Non-Basic and Basic variable
     * collections.
     * 
     * It returns a vector with the solutions for the initial non-basic variables.
     */
    public static double[] simplex(int[] N, int [] B, double[][]A, double[]b, double[]c){
        
        double v = initSimplex(N, B, A, b, c);
        
        //test to see if the solution is feasible, MIN_VALUE indicates infeasible
        if(v == Double.MIN_VALUE){
            System.out.println("The solution is infeasible");
            return null;
        }
        
        System.out.println();
        
        //the array of solutions for each X to return
        double []x = new double[c.length];
        
        v = simplexLoop(N, B, A, b, c, v);
        
        //if the solution is unbounded, -1 is returned
        if(v == Double.MAX_VALUE)
        {
            System.out.println("The solution is unbounded");
            return null;
        }
        
        //else an optimal solution is found
        //find the initial non-basic variables. These will have the numbers 1..N
        for(int i = 0; i<N.length; i++)
        {
            //for every initial non-basic variable, look to see if it is in the Basic variables
            for(int j = 0; j< B.length;j++)
            {
                //if this nonbasic variable is now a basic variable, set
                if(B[j] == i+1)
                {
                    //assign the b[j] value to the x[i] place in the return variable
                    //so the x_1 variable (now a basic variable) will be assigned to the 0 index
                    x[i] = b[j];
                }
            }
        }
        
        //assign objective function so that  it can be returned later
        objFn = c;
        
        //return the values for the objective function
        return x; 
    }
    /******************************************************************************************/
    /* 
     * simplexLoop is called by simplex and initSimplex to find the optimal solution to the 
     * problem. It is modularized here to eliminate code redundancy between the two methods. 
     * It takes the variable collections, Matrix, constraint and objective function vectors 
     * and a current solution to the maximization problem as arguments
     * 
     * It returns and updated solution to the maximization problem
     */
    private static double simplexLoop(int[] N, int [] B, double[][]A, double[]b, double[]c, double v)
    {
        
        int loopCounter = 0;
        //set entering variable to be "null" and init a finished variable
        int e = -1;
        boolean finished = false;
        
        //while some index j in N has cj > 0
        //choose a column to pivot on
        while(!finished)
        {
            loopCounter++;
            //reset the entering variable to be -1
            e = -1;
            
            //find a coeffiecient in the objective function that is positive
            for(int i = 0; i< c.length; i++)
            {
                if(c[i] > 0)
                {
                    e=i;
                    
                    //choose first positive coefficient
                    break; 
                }
            }
            
            //if we didn't find a coefficient that is positive, we are done
            if(e == -1)
            {
                finished = true;
            }
            
            if(!finished)
            {
                //find min ratio
                double minRatio = Double.MAX_VALUE;
                int l = -1;
                for(int i = 0; i< b.length; i++)
                {
                    //if the column value in matrix A is > 0
                    if(A[i][e] > 0){
                        
                        double delta = b[i]/A[i][e];

                        if(minRatio > delta)
                        {
                            minRatio = delta;
                            l = i;
                        }
                    }
                }
                
                //check for unbounded
                if(minRatio == Double.MAX_VALUE)
                {
                   return minRatio;
                }
                else{
                    v = pivot(N, B, A,b,c, v ,l ,e );
                }
            }
            
        }
        //System.out.println("loopCounter = " + loopCounter);
        return v;
    }
    /******************************************************************************************/
    /*
     * initSimplex takes the variable collections, Matrix, constraint and objective function vectors
     * and reports that the solution is infeasible, or it returns a basic feasible solution for
     * simplex to start with.
     * 
     * It returns the basic initial feasible solution.
     */
    private static double initSimplex(int[] N, int [] B, double[][]A, double[]b, double[]c){
        
        //get matrix dimensions
        int m = b.length; //this is the number of constraints
        int n = c.length; //this is the number of variables in the objective function
        
        //variable for the current solution to the objective function
        double v = 0;
        
        //init non basic vars, these are the variables in the objective function
        for(int i =0; i< n; i++)
        {
            N[i] = i + 1;
        }
        
        //init non basic vars, these are the slack variables
        for(int i =0; i< m; i++)
        {
            //set basic variable subscript to n + i + 1, n + 1 ... n + m
            B[i] = n + i + 1;
        }
        
        //find the smalles b value index
        int k = findMinb(b);
        
        //do we have a feasible solution with no negative b? if so, return to simplex
       if(b[k] >= 0)return v;
        
        //otherwise construct the aux function by adding -x_0 to each LHS and set the objective function to -x_0
        int l = k; //the basic variable to leave in the pivot operation

        //re-init the Non basic vars to include -x_0
        int[] N_ = new int[n+1];
        for(int i =0; i< n + 1; i++)
        {
            N_[i] = i;
        }
        
        //init new objective function to -x_0, all other vars have a coef of 0
        double[] c_ = new double[c.length + 1];
        c_[0]=-1;
        
        //add x_0 to the matrix
        double A_[][] = new double[A.length][A[0].length + 1];
        
        for(int i = 0; i< A.length; i++)
        {
            //add the x_0 term at the index 0
            A_[i][0] = -1;
            
            //fill in the rest of the matrix, one place higher
            for(int j = 0; j< A[0].length; j++)
            {
                A_[i][j+1] = A[i][j];
            }
        }
        
        //pivot the solution
        v = pivot(N_, B, A_, b, c_, v, l, 0);
        
        //iterate simplex
        v = simplexLoop(N_, B, A_, b, c_, v);
        
        //if the optimal solutions for -x_0 is 0
        if(v==0)
        {
            //see if x_0 is basic
            for(int i = 0; i< B.length; i++)
            {
                if(B[i] == 0)
                {
                    int pivotIdx = -1;
                    //perform pivot to make it non basic
                    //pivot with any j in N, so that a_[i][j] != 0
                    for(int j = 0; j< N.length; j++)
                    {
                        if(A_[i][j] != 0)
                        {
                          //first non-zero non-basic variable index
                          pivotIdx = j;
                          break;
                        }
                    }
                    v = pivot(N_, B, A_, b, c_,v, i, pivotIdx);
                }
            }
            
            //remove x_0 from the non basic vars and constraints      
            int nIdx = 0;
            for(int i = 0; i< N_.length; i++)//number of vars or columns
            {
                //this is not x0
                if(N_[i] != 0)
                {
                    //add the variable to N
                    N[nIdx] = N_[i];

                    //set the new matrix value
                    for(int j=0; j< B.length;j++) //number of constraints or rows
                    {
                        A[j][nIdx] = A_[j][i];
                    }

                    //increase the index for N when you add a variable to it
                    //this will match i until x0 is encountered, the it will be i-1
                    nIdx++;
                }
            }
            
            //new matrix to hold the new C coefficients
            double[][] tMatrix = new double[B.length][N.length];
            
            //replace each basic var in obj function by RHS of associated constraint
            for(int i = 0; i<B.length; i++)
            {
                //we know original vars were 1..N, so look for those in B
                if(B[i] <= n)
                {
                    int bIdx = B[i]-1; //vars are 1..N and indexes are 0..N-1                
                    
                    //multiply coef from original NB var times the new constraint coef.
                    //i is the row, j is the column
                    for(int j = 0; j< N.length; j++)
                    {
                        //new matrix for coeficient multiplication
                        tMatrix[i][j] = (c[bIdx] * A[i][j])*-1; //add negative one because these vars are on the RHS
                        
                    }
                  
                    v += (c[bIdx]*b[i]);
                }
            }
            
            //initialize a new vector to hold the new calculations for the c coefficients
            c_ = new double[N.length];
            
            //after computing the matrix, combine like terms to find new coef for C
            for(int i=0; i< B.length; i++)
            {
                for(int j = 0; j< N.length; j++)
                {
                    c_[j] += tMatrix[i][j];

                }
            }
         
            //if j is in the non basic variables, combine like terms
            //count from 0..N
            for(int j = 0; j< N.length; j++)
            {
                //loop through the Non Basic vars
                for(int h=0; h<N.length; h++)
                {
                    //if this basic bar is a non basic var in the original objective function
                    if(N[h] == j+1)
                    {
                        //c[j] is the variable's original position in the objective function
                        //add the value from the original objective function to the new objective function
                        //combining like terms
                        c_[h] += c[j];
                    }
             
                }
            }
            //return new modified slack form of L 
            System.arraycopy(c_, 0, c, 0, c.length);
        }
        
        else
        {
           //return min value to indicate the solution is infeasible
           return Double.MIN_VALUE;
        }
        return v;
    }
    /******************************************************************************************/
    /*
     * pivot takes the variable collections, Matrix, constraint and objective function vectors, 
     * entering and leaving variable, and a current solution.
     * 
     * It returns a pivoted solution for simplex to process working towards optimality.
     */
    private static double pivot(int[] N, int [] B, double[][]A, double[]b, double[]c, double v, int l, int e)
    {
        //length of the nonbasic variables
        int n = N.length;
        
        //length of the basic variable
        int m = B.length;
        
        //create new matrix
        double[][]A_ = new double[A.length][A[0].length];
        
        //create new b vector the length of the Basic variables
        double[]b_ = new double[B.length];
        
        //create new c vector the length of the Non-basic variables
        //double[]c_ = new double[c.length];
        double[]c_ = new double[N.length];
        
        //create new N and B
        int[]N_ = new int[N.length];
        int[]B_ = new int[B.length];
        
        //create variable for new v
        double v_;
        
        //compute the new b for the leaving variable
        b_[l] = b[l]/A[l][e];
        
        //for the new leaving variable, compute the coefficients
        for(int i = 0; i< N.length; i++)
        {
            //ignore the entering variable
            if(i != e)
            {
                //compute coefficient ratio 
                A_[l][i] = A[l][i]/A[l][e];
            }
            //if this is the entering var, set the coef to it's inverse as the slack value is 1 
            //and we are dividing by the value of the entering coefficient
            else
            {
                A_[l][i] = 1/A[l][i];
            }
        }  
        
        //compute coefficients for the remaining constraints
        for(int i =0; i< b.length; i++)
        {
            //ignore the leaving variable row we calculated above
            if(i != l)
            {
                //the new b is the old b for this row - the leaving b * the number in the entering vector for this row
                b_[i] = b[i] - (A[i][e] * b_[l]);
                
                 for(int j = 0; j< N.length; j++)
                {
                    //ignore the entering vector variable, this will be 0
                    if(j != e)
                    {
                        //new vector value calculation
                       //the current value - (the vector value for this row * new leaving value for this column vector )
                       A_[i][j]=A[i][j] - (A[i][e]*A_[l][j]);
                    }
                    //new vector calculation
                    //0 - the exisiing value * the new leaving value for this column vector
                    else
                    {
                        A_[i][e]= 0 - (A[i][e]*A_[l][j]);
                    }
                }
            }
        }
        
        //in the book, the algorithm specifies the new b of the entering var, but I think the leaving var is meant
        //compute the new optimal value for z
        v_ = v + (c[e]*b_[l]);
        
        //compute the obective function coefficients
        for(int j = 0; j< N.length; j++)
        {
            if(j != e)
            {
                //the coefficient for the staying var is itself, minus 
                //the coefficient of the entering variable *  the new value for this column in the leaving var
                c_[j] = c[j]-(c[e]*A_[l][j]);
                
            }
            //the coefficient for the leaving variable is going to be 0 - new value of the leaving vector 
            // * the current coefficient for the entering vector
            else
            {
                c_[j]= 0 - (c[e]*A_[l][e]);
            }
        }
        
        
        //compute new B and N vectors
        for(int i = 0; i< B.length; i++)
        {
            //if it isn't leaving,the B stays the same
            if(i != l)
            {
                B_[i] = B[i];
            }
            //if it is, then assign the entering value in its place
            else
            {
                 B_[i]= N[e];
            }
        }
        
        for(int i = 0; i< N.length; i++)
        {
            //if it isn't entering,the N stays the same
             if(i != e)
            {
                N_[i] = N[i];
            }
            else
            {
                //if it is, then assign the leaving value in its place
                N_[i] = B[l];
            }
            
        }
        //assign returning variables
        System.arraycopy(N_, 0, N, 0, N.length);
        System.arraycopy(B_, 0, B, 0, B.length);
        System.arraycopy(A_, 0, A, 0, A.length);
        System.arraycopy(b_, 0, b, 0, b.length);
        System.arraycopy(c_,0, c, 0, c.length);
       
        return v_;

    }
    /******************************************************************************************/
    /*
     * print is a helper function for debugging code. It prints out the current status of the solution 
     * to the LP to the standard output. 
     */
    static void print(int[] N, int [] B, double[][]A, double[]b, double[]c, double v, int l, int e)
    {
        System.out.print("N: ");
        for(int i = 0; i< N.length; i++)
        {
            System.out.print(N[i] + ", ");
        }
        System.out.println();
        
        System.out.print("B: ");
        for(int i = 0; i< B.length; i++)
        {
            System.out.print(B[i] + ", ");
        }
        System.out.println();
        
        System.out.print("A: ");
        for(int i = 0; i< A.length; i++)
        {
            for(int j = 0; j<A[0].length; j++)
            {
                System.out.print(A[i][j] + " ");
            }
              System.out.println();
        }
        
        System.out.print("b: ");
        for(int i = 0; i< b.length; i++)
        {
            System.out.print(b[i] + ", ");
        }
        System.out.println();
        
        System.out.print("c: ");
        for(int i = 0; i< c.length; i++)
        {
            System.out.print(c[i] + ", ");
        }
        System.out.println();
        
        System.out.println("v = " + v );
        System.out.println("l= " + l );
        System.out.println("e= " + e);
    }
    /******************************************************************************************/
    /*
     * findMinb is a method that is used to determine which constraint has the least value
     * It is used in the initSimplex method to determine if the basic solution is feasible
     * 
     * It returns the index of the minimum b value.
     */
    private static int findMinb(double[]b){
        
        //set a variable to hold the min b value and min index
        double min = Double.MAX_VALUE;
        int minIdx = -1;
        
        //loop through all the b values, select the min value
        for(int i=0; i< b.length; i++){
            if(b[i] < min){
                min=b[i];
                minIdx = i;
            }
        }
        
        return minIdx;
    }
    /******************************************************************************************/
}
