package sph9_iga;

import static java.lang.Math.sqrt;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;

public class SPH extends Thread{
	//constants
	private final int FRAME=1000;
	private final double SCALE=3;
	private final int RS=1;
	
	//fields
	private Parameter p=new Parameter();
	private Particle[] part=new Particle[p.POP];
	
	//density fields
	private int xm=(int)((p.MAX[0]-p.MIN[0])*SCALE);
	private int ym=(int)((p.MAX[1]-p.MIN[1])*SCALE);
	private int zm=(int)((p.MAX[2]-p.MIN[2])*SCALE);
	private double[][][] dens=new double[zm][zm][zm];
	private byte[][][] densB=new byte[zm][zm][zm];
	
	private int xd=(int)((p.MAX[0]-p.MIN[0])*SCALE/(p.H*RS/p.SCALE))+1;
	private int yd=(int)((p.MAX[1]-p.MIN[1])*SCALE/(p.H*RS/p.SCALE))+1;
	private int zd=(int)((p.MAX[2]-p.MIN[2])*SCALE/(p.H*RS/p.SCALE))+1;
	private Room[][][] room=new Room[xd][yd][zd];
	private int members=10;
	
	//constructor
	SPH(){
		//initialize room
		for(int i=0;i<xd;i++){
			for(int j=0;j<yd;j++){
				for(int k=0;k<zd;k++){
					room[i][j][k]=new Room(members);
				}
			}
		}
	}
	
	//methods
	
	//initialize
	private void init(){
		for(int i=0;i<p.POP;i++){
			part[i]=new Particle(p.DIM,p.NB);
		}
		
		int count=0;
		for(int x=0;x<p.PX[0];x++){
			for(int y=0;y<p.PX[1];y++){
				for(int z=0;z<p.PX[2];z++){
					part[count].loc[0]=(x+0.5)*p.D;
					part[count].loc[1]=(y+0.5)*p.D;
					part[count].loc[2]=(z+0.5)*p.D;
					count++;
				}
			}
		}
		
		for(int x=0;x<p.PX[0];x++){
			for(int y=0;y<p.PX[1];y++){
				for(int z=0;z<p.PX[2];z++){
					part[count].loc[0]=60-(x+0.5)*p.D;
					part[count].loc[1]=10-(y+0.5)*p.D;
					part[count].loc[2]=(z+0.5)*p.D;
					count++;
				}
			}
		}
		
		createFile(0);
	}
	
	//calculate density and pressure
	private void calcDP(){
		double r2,sum,dummy;
		int num=0;
		
		for(int i=0;i<p.POP;i++){
			sum=0;
			
			for(int j=0;j<part[i].nbCount;j++){
				num=part[i].nbList[j];
				r2=norm2(part[i].loc,part[num].loc)*p.SCALE*p.SCALE;
				
				if(r2<p.H2){
					dummy=p.H2-r2;
					sum+=dummy*dummy*dummy;
				}
			}//for j
			
			sum+=p.H2*p.H2*p.H2;
			
			part[i].dens=p.PMASS*p.POLY6KERN*sum;
			part[i].pres=p.I_STIFF*(part[i].dens-p.R_DENS);
		}//for i
		
	}
	
	//calculate acceleration and color field
	private void calcForce(){
		double pterm,vterm,r,color,dummy,dummy2;
		double[] pacc=new double[p.DIM];
		double[] vacc=new double[p.DIM];
		int num=0;
		
		//debug variable
		double min=0;
		double max=0;
				
		for(int i=0;i<p.POP;i++){
			for(int d=0;d<p.DIM;d++){
				pacc[d]=0;
				vacc[d]=0;
			}
			color=0;
			
			for(int j=0;j<part[i].nbCount;j++){
				num=part[i].nbList[j];
				
				if(i!=num){
					r=sqrt(norm2(part[i].loc,part[num].loc))*p.SCALE;
					
					if(r<p.H){
						dummy=p.H-r;
						pterm=dummy*dummy*(part[i].pres+part[num].pres)/(r*part[num].dens);
						vterm=dummy/part[num].dens;
						dummy2=p.H2-r*r;
						color+=1.0/part[num].dens*dummy2*dummy2*dummy2;
						for(int d=0;d<p.DIM;d++){
							pacc[d]+=pterm*(part[i].loc[d]-part[num].loc[d]);
							vacc[d]+=vterm*(part[num].vel[d]-part[i].vel[d]);
						}
					}
				}//if j
			}//for j
			
			part[i].color=p.PMASS*p.POLY6KERN*color;
			for(int d=0;d<p.DIM;d++){
				part[i].acc[d]=p.PMASS/part[i].dens*(-0.5*p.SPIKY_KERN*pacc[d]
						+p.VISCOSITY_KERN*p.VISC*vacc[d]);
			}
			
			//debug
			if(i==0){
				min=max=part[i].color;
			}else{
				dummy=part[i].color;
				if(dummy<min){
					min=dummy;
				}
				if(max<dummy){
					max=dummy;
				}
			}
		}//for i	
		
		System.out.println("color max:"+max+",min:"+min);
	}
	
	//calculate location
	private void calcLoc(){
		double dummy;
		
		for(int i=0;i<p.POP;i++){
			//acceleration limit
			dummy=norm2(part[i].acc,new double[]{0,0,0});
			if(p.ALIM*p.ALIM<dummy){
				dummy=sqrt(dummy);
				for(int d=0;d<p.DIM;d++){
					part[i].acc[d]*=p.ALIM/dummy;
				}
			}
			
			for(int d=0;d<p.DIM;d++){
				//lower boundary condition
				dummy=2*p.RAD-(part[i].loc[d]-p.MIN[d])*p.SCALE;
				if(0.000001<dummy){
					part[i].acc[d]+=p.E_STIFF*dummy-p.E_DAMP*part[i].vel[d];
				}
				
				//upper boundary condition
				dummy=2*p.RAD-(p.MAX[d]-part[i].loc[d])*p.SCALE;
				if(0.000001<dummy){
					part[i].acc[d]+=-p.E_STIFF*dummy-p.E_DAMP*part[i].vel[d];
				}
			}
			
			part[i].acc[2]+=p.gravity;
			
			for(int d=0;d<p.DIM;d++){
				part[i].vel[d]+=part[i].acc[d]*p.DT;
				part[i].loc[d]+=part[i].vel[d]*p.DT/p.SCALE;
			}
		}//for i
	}
	
	//sub calculation methods
	private double norm2(double[] x,double[] y){
		double dummy=0;
		for(int d=0;d<x.length;d++){
			dummy+=(x[d]-y[d])*(x[d]-y[d]);
		}
		return dummy;
	}
	
	//check neighborhood
	private void checkNb(){
		double dr2=0;
		double lh2=4*p.H2;
		
		for(int i=0;i<p.POP;i++){
			part[i].nbCount=0;
			for(int j=0;j<p.NB;j++){
				part[i].nbList[j]=0;
			}
		}
		
		for(int i=0;i<p.POP;i++){
			for(int j=(i+1);j<p.POP;j++){				
				if(part[i].nbCount<p.NB){
					//calculate distance
					dr2=0;
					for(int d=0;d<p.DIM;d++){
						dr2+=(part[i].loc[d]-part[j].loc[d])*(part[i].loc[d]-part[j].loc[d]);
					}
					dr2*=p.SCALE*p.SCALE;
					
					//list neighborhood particle
					if(dr2<lh2){
						part[i].nbList[part[i].nbCount]=j;
						part[i].nbCount++;
						if(part[j].nbCount<p.NB && i!=j){
							part[j].nbList[part[j].nbCount]=i;
							part[j].nbCount++;	
						}
					}
				}//if count
			}//for j
		}//for i
	}
	
	//create directory and location file
	private void createDirectory(){
		//create directory
		File dir=new File(System.getProperty("user.dir")+"\\SPH"+p.NUM);
		if(!dir.exists()){
			dir.mkdir();	
		}
		
		//create density directory
		File dirD=new File(System.getProperty("user.dir")+"\\SPH_DF3"+p.NUM);
		if(!dirD.exists()){
			dirD.mkdir();	
		}
	}
	
	private void createFile(int num){
		//create data file
		try{
			File dataFile=new File("SPH"+p.NUM+"\\sph"+num+".dat");
			if(dataFile.exists()){
				dataFile.delete();
			}
			dataFile.createNewFile();
			FileWriter fw=new FileWriter(dataFile,true);
			
			for(int i=0;i<p.POP;i++){
				fw.write("<"+(float)part[i].loc[0]+","
							+(float)part[i].loc[2]+","
							+(float)part[i].loc[1]+">, ");
			}
			
			fw.close();
			
			System.out.println("Frame:"+num);
			System.out.println("dens:"+part[0].dens+",pres:"+part[0].pres+",nb:"+part[0].nbCount
					+",x:"+part[0].loc[2]+",v:"+part[0].vel[2]+",a:"+part[0].acc[2]);
			
			//density file
			createDF3(num);
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	//create df3 file
	private void createDF3(int num){
		try{
			//calculate density file
			double min,max,sum;
			min=max=0;
			int x,y,z;
			
			//initialize member of each room
			for(int i=0;i<xd;i++){
				for(int j=0;j<yd;j++){
					for(int k=0;k<zd;k++){
						room[i][j][k].num=0;
					}
				}
			}
			
			//count member
			for(int i=0;i<p.POP;i++){
				x=(int)(part[i].loc[0]*SCALE/(p.H*RS/p.SCALE));
				y=(int)(part[i].loc[1]*SCALE/(p.H*RS/p.SCALE));
				z=(int)(part[i].loc[2]*SCALE/(p.H*RS/p.SCALE));
				//problem
				if(-1<x && -1<y && -1<z && x<xd && y<yd && z<zd){
					if(room[x][y][z].num<members){
						room[x][y][z].member[room[x][y][z].num]=i;
						room[x][y][z].num++;
					}else{
						System.out.println("out");
					}
				}
			}
			
			//calculate density at each point
			for(int i=0;i<xm;i++){
				for(int j=0;j<ym;j++){
					for(int k=0;k<zm;k++){
						sum=0;
						x=(int)(i/(p.H*RS/p.SCALE));
						y=(int)(j/(p.H*RS/p.SCALE));
						z=(int)(k/(p.H*RS/p.SCALE));
						
						//origin
						sum+=checkSum(new int[]{x,y,z},new int[]{i,j,k});
						
						//1 basis
						if((x+1)<xd){
							sum+=checkSum(new int[]{x+1,y,z},new int[]{i,j,k});
						}
						if(0<(x-1)){
							sum+=checkSum(new int[]{x-1,y,z},new int[]{i,j,k});
						}
						if((y+1)<yd){
							sum+=checkSum(new int[]{x,y+1,z},new int[]{i,j,k});
						}
						if(0<(y-1)){
							sum+=checkSum(new int[]{x,y-1,z},new int[]{i,j,k});
						}
						if((z+1)<xd){
							sum+=checkSum(new int[]{x,y,z+1},new int[]{i,j,k});
						}
						if(0<(z-1)){
							sum+=checkSum(new int[]{x,y,z-1},new int[]{i,j,k});
						}
						
						//2basis
						if((x+1)<xd && (y+1)<yd){
							sum+=checkSum(new int[]{x+1,y+1,z},new int[]{i,j,k});
						}
						if((y+1)<yd && (z+1)<zd){
							sum+=checkSum(new int[]{x,y+1,z+1},new int[]{i,j,k});
						}
						if((z+1)<zd && (x+1)<xd){
							sum+=checkSum(new int[]{x+1,y,z+1},new int[]{i,j,k});
						}
						
						if((x+1)<xd && 0<(y-1)){
							sum+=checkSum(new int[]{x+1,y-1,z},new int[]{i,j,k});
						}
						if((y+1)<yd && 0<(z-1)){
							sum+=checkSum(new int[]{x,y+1,z-1},new int[]{i,j,k});
						}
						if((z+1)<zd && 0<(x-1)){
							sum+=checkSum(new int[]{x-1,y,z+1},new int[]{i,j,k});
						}
						
						if(0<(x-1) && (y+1)<yd){
							sum+=checkSum(new int[]{x-1,y+1,z},new int[]{i,j,k});
						}
						if(0<(y-1) && (z+1)<zd){
							sum+=checkSum(new int[]{x,y-1,z+1},new int[]{i,j,k});
						}
						if(0<(z-1) && (x+1)<xd){
							sum+=checkSum(new int[]{x+1,y,z-1},new int[]{i,j,k});
						}
						
						if(0<(x-1) && 0<(y-1)){
							sum+=checkSum(new int[]{x-1,y-1,z},new int[]{i,j,k});
						}
						if(0<(y-1) && 0<(z-1)){
							sum+=checkSum(new int[]{x,y-1,z-1},new int[]{i,j,k});
						}
						if(0<(z-1) && 0<(x-1)){
							sum+=checkSum(new int[]{x-1,y,z-1},new int[]{i,j,k});
						}
						
						//3 basis
						if((x+1)<xd && (y+1)<yd && (z+1)<zd){
							sum+=checkSum(new int[]{x+1,y+1,z+1},new int[]{i,j,k});
						}
						if(0<(x-1) && (y+1)<yd && (z+1)<zd){
							sum+=checkSum(new int[]{x-1,y+1,z+1},new int[]{i,j,k});
						}
						if((x+1)<xd && 0<(y-1) && (z+1)<zd){
							sum+=checkSum(new int[]{x+1,y-1,z+1},new int[]{i,j,k});
						}
						if((x+1)<xd && (y+1)<yd && 0<(z-1)){
							sum+=checkSum(new int[]{x+1,y+1,z-1},new int[]{i,j,k});
						}
						if((x+1)<xd && 0<(y-1) && 0<(z-1)){
							sum+=checkSum(new int[]{x+1,y-1,z-1},new int[]{i,j,k});
						}
						if(0<(x-1) && (y+1)<yd && 0<(z-1)){
							sum+=checkSum(new int[]{x-1,y+1,z-1},new int[]{i,j,k});
						}
						if(0<(x-1) && 0<(y-1) && (z+1)<zd){
							sum+=checkSum(new int[]{x-1,y-1,z+1},new int[]{i,j,k});
						}
						if(0<(x-1) && 0<(y-1) && 0<(z-1)){
							sum+=checkSum(new int[]{x-1,y-1,z-1},new int[]{i,j,k});
						}
						
						dens[i][j][k]=sum;
						
						//check max
						if(i==0 && j==0 && k==0){
							min=max=sum;
						}else{
							if(sum<min){
								min=sum;
							}
							if(max<sum){
								max=sum;
							}
						}//if max or min
					}
				}
			}
			
			System.out.println("dens max:"+max+",min:"+min);
			
			//create data file
			File dataFile=new File("SPH_DF3"+p.NUM+"\\sph"+num+".df3");
			if(dataFile.exists()){
				dataFile.delete();
			}
			dataFile.createNewFile();
			FileOutputStream fos=new FileOutputStream(dataFile);
			
			//write DF3 file
			
			//header
			short[] size=new short[]{(short)zm,(short)zm,(short)zm};
			byte dummy1,dummy2;
			for(int d=0;d<size.length;d++){
				dummy1=(byte)(size[d]>>8);
				dummy2=(byte)(size[d]);
				fos.write(dummy1);
				fos.write(dummy2);
			}
			
			for(int i=0;i<zm;i++){
				for(int j=0;j<zm;j++){
					for(int k=0;k<zm;k++){
						densB[i][j][k]=(byte)(-128*(dens[k][j][i]-min)/(max-min)+1);
					}
				}
			}
			
			for(int i=0;i<zm;i++){
				for(int j=0;j<zm;j++){
					fos.write(densB[i][j]);
				}
			}
			
			fos.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	//check sum of density
	private double checkSum(int[] x,int[] xp){
		double sum=0;
		double dr2,dummy;
		int num;
		
		for(int i=0;i<room[x[0]][x[1]][x[2]].num;i++){
			dr2=0;
			num=room[x[0]][x[1]][x[2]].member[i];
			for(int d=0;d<p.DIM;d++){
				dr2+=(xp[d]/SCALE-part[num].loc[d])*(xp[d]/SCALE-part[num].loc[d]);
			}
			dr2*=p.SCALE*p.SCALE;
			
			if(dr2<p.H2){
				dummy=p.H2-dr2;
				sum+=dummy*dummy*dummy;
			}
		}
		
		return sum*100000;
	}
		
	//run method
	public void run(){
		createDirectory();
		init();
		checkNb();
		for(int i=1;i<FRAME;i++){
			if(i%5==0){
				checkNb();
			}
			calcDP();
			calcForce();
			calcLoc();
			createFile(i);
		}
	}
}
