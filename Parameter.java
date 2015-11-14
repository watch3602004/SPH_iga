package sph9_iga;

import static java.lang.Math.*;

public class Parameter {
	//constants
	protected final int NUM=0;
	
	//simulation parameters
	protected final int DIM=3;
	protected final double SCALE=0.04;
	
	//external conditions
	protected final double[] INIMIN=new double[]{0,0,0};
	protected final double[] INIMAX=new double[]{10,10,10};
	protected final double[] MIN=new double[]{0,0,0};
	protected final double[] MAX=new double[]{60,10,60};
	protected final double gravity=-9.8;
	
	//SPH parameters
	protected final double PMASS=0.010543;//kg 0.020543
	protected final double R_DENS=1000;//kg.m-3
	protected final double PDIST=pow(PMASS/R_DENS,1/3.0);//m
	protected final double D=PDIST*0.9/SCALE;//why 0.87?
	
	protected final int[] PX=new int[DIM];
	protected final int POP;
	
	protected final double H=0.03;//m 0.1
	protected final double H2=H*H;
	
	protected final double I_STIFF=1;
	protected final double E_STIFF=10000;
	protected final double VISC=0.05;
	protected final double ALIM=200;//m.s-2
	protected final double RAD=0.004;//m
	protected final double E_DAMP=256;
	
	protected final double DT=0.004;
	
	//smoothing kernels
	protected final double POLY6KERN=315/(64*PI*pow(H,9));
	protected final double SPIKY_KERN=-45/(PI*pow(H,6));
	protected final double VISCOSITY_KERN=45/(PI*pow(H,6));
	
	//neighbor list
	protected final int NB=200;
	
	//fields

	//constructor
	Parameter(){
		int pop=1;
		
		for(int d=0;d<DIM;d++){
			PX[d]=(int)((INIMAX[d]-INIMIN[d])/D)+1;
			pop*=PX[d];
		}
		
		POP=2*pop;
		System.out.println("POP:"+POP+",H/PDIST:"+H/PDIST+",D"+D);
	}
	
	//methods
}
