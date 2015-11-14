package sph9_iga;

public class Particle {
	//constants
	
	//fields
	protected double[] loc;
	protected double[] vel;
	protected double[] acc;
	protected double[] gclv;//gradient of smoothed color field

	protected double dens;
	protected double pres;
	protected double color;//smoothed color field
	protected double gcl;//absolute of gradient of smoothed color field
	protected double lcl;//Laplacian of smoothed color field
	
	protected int[] nbList;
	protected int nbCount=0;
	
	//constructor
	Particle(int DIM,int NB){
		loc=new double[DIM];
		vel=new double[DIM];
		acc=new double[DIM];
		gclv=new double[DIM];
		nbList=new int[NB];
	}
	
	//methods
	
	//initialize
	public void init(){
		for(int i=0;i<loc.length;i++){
			loc[i]=0;
			vel[i]=0;
			acc[i]=0;
		}
	}
}
