
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Random;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author shannon
 */
public class CreateLPPInput {
    
    public static void main(String args[])
    {
     try
        {
            FileWriter writer = new FileWriter(args[0]);
            BufferedWriter buffer = new BufferedWriter(writer);
            
            Random r = new Random(System.currentTimeMillis());
        
         
            int size1 = Integer.parseInt(args[1]);
            int size2 = Integer.parseInt(args[2]);
            
            //loop through all but the 0 and last, which are our added source and sink
            
            for(int i = 0; i < size1; i++){
				buffer.write(r.nextInt(50) + " "); 
			}
			buffer.write("\n");   
			
			for(int j = 0; j < size2; j++){
				buffer.write(r.nextInt(50) + " "); 
			}
			buffer.write("\n");   
            
            for(int i =0; i< size2;i++)
            {
                //buffer.write(0 + " ");
                
                for(int j =0; j< size1;j++)
                {
                    buffer.write(r.nextInt(50) + " ");   
                }
                buffer.write("\n");   
            }
            
            //write last row for sink
             //for(int j =1; j< size+1;j++)
            //{
               //buffer.write(0 + " ");
           // }
             
            buffer.close();
        }
        catch(Exception e)
        {
            System.out.println("An error occurred " + e.getMessage());
            
        }
    }

}
