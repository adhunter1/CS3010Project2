import java.io.*;
import java.util.*;

public class assignment2 {

	public static void main(String[] args) {
		int num=0;
		double[][] coeff=new double[1][1];
		double[] sol= new double[1];
		String[] ff=new String [3];
		Scanner scan = new Scanner(System.in);
		do {
			System.out.println("Please enter a valid file name. Ex. \"Gaussian sys1.lin\" or \"Gaussian --spp sys1.lin\"");
			String fileName=scan.nextLine();
			ff =fileName.split(" ");
		}while(!ff[0].equalsIgnoreCase("Gaussian"));
		
		
		
		try {
			String fileName;
			if(ff[1].equalsIgnoreCase("--spp")){
				fileName=ff[2];

			}else {
				fileName=ff[1];
			}

			File file = new File("src//assignment2CS3010//"+fileName);
					//"C:\\Users\\Alex\\eclipse-workspace\\assignment2CS3010\\src\\assignment2CS3010\\sys1.lin");
			BufferedReader br = new BufferedReader(new FileReader(file));
			StringBuffer sB = new StringBuffer();
			String line;
			//set num to first number
			num=Integer.parseInt(br.readLine());
			coeff=new double[num][num];
			sol=new double[num];
			//go thru a loop till num
			for(int i= 0; i<num;i++)
			{
				line=br.readLine();
				String[] eachLine=line.split("\\s+");
				//reads in each line then splits in and places it into array
				for(int j=0;j<num;j++) {
					coeff[i][j]=Double.parseDouble(eachLine[j]);
				}
			}
			line=br.readLine();
			
			String[] eachLine=line.split("\\s+");
			for(int i=0;i<num;i++) {
				sol[i]=Double.parseDouble(eachLine[i]);
			}
			
			

			br.close();
			System.out.println(sB.toString());
		} catch (IOException e) {
			e.printStackTrace();
		}
		if(ff[1].equalsIgnoreCase("--spp")){
			//then do the scaled partial pivoting
			sol=SSPGaussian(coeff,sol);
		}else {
			//do reg
			sol=NaiveGaussianElimination(coeff,sol);
		}
		
		//write in sys1.sol
		String fileWrite="src\\assignment2CS3010\\sys1.sol";
		BufferedWriter bw= null;
		FileWriter fw = null;
		try {
			fw=new FileWriter(fileWrite);
			bw= new BufferedWriter(fw);
			for(int i =0;i<num;i++) {
				bw.write(sol[i]+" ");
			}
		}catch(IOException e) {
			e.printStackTrace();
		}finally {
			try {

				if (bw != null)
					bw.close();

				if (fw != null)
					fw.close();

			} catch (IOException e) {

				e.printStackTrace();

			}
		}
		
	}
	
	public static double[] FwdElimination(double[][] coeff, double[] constant) {
		int n = coeff.length;
		for(int i = 0 ; i < n; i++){
            for (int k = i + 1; k < n; k++){
                double scale = -(coeff[k][i])/(coeff[i][i]);

                for(int j = i; j < n; j++){
                    if (i == j){
                        coeff[k][j] = 0;
                    }
                    else{
                        coeff[k][j] += scale * coeff[i][j];
                    }
                }

                constant[k] += scale * constant[i];
            }
        }
		return constant;
	}
	
	public static double[] BackSubstitution(double[][] coeff, double[] constant, double [] sol){

		 

        for (int i = coeff.length - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < coeff.length; j++) {
                sum += coeff[i][j] * sol[j];
            }
            sol[i] = (constant[i] - sum) / coeff[i][i];
        }


        return sol;
    }

	 public static double[] NaiveGaussianElimination(double[][] coeff, double[] constant){

	        int n = coeff.length;
	        double [] sol=new double[n];
	        sol=FwdElimination(coeff,constant);
	        return BackSubstitution(coeff, constant,sol);
	    }

	    //SCALED PARTIAL PIVOTING.
	 
	 public static double[] SPPFwdElimination(double[][] coeff, double[] constant, int[] ind ){
	     int n=coeff.length;
		 double[] scaling= new double [n];
		 for(int i=0;i<n;i++) {
			 double smax=0;
			 for(int j=0;j<n;j++)
			 {
				 smax=Math.max(smax,Math.abs(coeff[i][j]));
			 }
			 scaling[i]=smax;
		 }
		 for(int k=0;k<n-1;k++) {
			 double rmax=0;
			 int maxInd=k;
			 for(int i=k;i<n;i++)
			 {
				 double r=Math.abs(coeff[(ind[i])][k]/scaling[ind[i]]);
				 if(r>rmax) {
					 rmax=r;
					 maxInd=i;
				 }
			 }
			 int temp=ind[maxInd];
			 ind[maxInd]=ind[k];
			 ind[k]=temp;
			 for(int i=k+1;i<n;i++) {
				 double mult=coeff[ind[i]][k]/coeff[ind[k]][k];
				 for(int j=k+1;j<n;j++) {
					 coeff[ind[i]][j]=coeff[ind[i]][j]-mult*coeff[ind[k]][j];
				 }
				 constant[ind[i]]=constant[ind[i]]-mult*constant[ind[k]];
			 }
			 
		 }
		 return constant;
	 }
	 public static double[] SPPBackSubst(double[][] coeff, double[] constant, double[] sol, int[] order){
	        double[] answer = new double[coeff.length];
	        int last = order[order.length - 1];


	        answer[answer.length - 1] = constant[last]/coeff[last][coeff[last].length - 1];

	        for(int i = order.length - 2; i >= 0; i--){
	            double sum = 0;
	            int currentRow = order[i];

	            for(int j = i + 1; j < coeff[i].length; j++){
	                sum += coeff[currentRow][j] * answer[j];
	            }

	            answer[i] = (constant[currentRow] - sum)/ coeff[currentRow][i];
	        }

	        return answer;
	    }
	 public static double[] SSPGaussian(double[][] coeff, double[] constant) {
		 int n=coeff.length;
		 double[] sol= new double[n];
		 int[] ind= new int[n];
		 for(int i=0;i<n;i++) {
			 ind[i]=i;
		 }
		 sol=SPPFwdElimination(coeff,constant,ind);
		 sol=SPPBackSubst(coeff,constant,sol,ind);
		 
		 return sol;
	 }
