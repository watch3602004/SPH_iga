package sph9_iga;

import java.io.File;
import java.io.FileOutputStream;

public class DF3 {
	//constant
	private final int NUM;
	
	//fields
	
	//constructor
	DF3(int num){
		this.NUM=num;
		
		//create directory
		File dir=new File(System.getProperty("user.dir")+"\\SPH_DF3"+num);
		if(!dir.exists()){
			dir.mkdir();	
		}
	}
	
	//methods
	public void createDF3(int x,int y,int z,byte[][][] dens,int num){
		try{
			//create data file
			File dataFile=new File("SPH_DF3"+NUM+"\\sph"+num+".df3");
			if(dataFile.exists()){
				dataFile.delete();
			}
			dataFile.createNewFile();
			FileOutputStream fos=new FileOutputStream(dataFile);
			
			//write DF3 file
			
			//header
			short[] size=new short[]{(short)x,(short)y,(short)z};
			byte dummy1,dummy2;
			for(int d=0;d<size.length;d++){
				dummy1=(byte)(size[d]>>8);
				dummy2=(byte)(size[d]);
				fos.write(dummy1);
				fos.write(dummy2);
			}
			
			//density data
			for(int i=0;i<x;i++){
				for(int j=0;j<y;j++){
					for(int k=0;k<z;k++){
						fos.write(dens[i][j][k]);
					}
				}
			}
			
			fos.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
}
