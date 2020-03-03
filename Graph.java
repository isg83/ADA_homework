import java.util.*;
import java.io.*;

public class Graph{
//metodo devuelve objeto HashMap con N nodos, argumento N
  public static HashMap <String,Double> node(int N){
    //creacion de objeto de la clase HashMap llamado N_nodes
    HashMap <String,Double> N_nodes = new HashMap <String,Double>();
    //variable de metodo para contener el número de Nodos
    int nodes = N;
    double x=0;
    //ciclo para meter nodo a la lista HashMap
    for(int i=0;i< nodes; i++) {
	x=Math.random();
      //si el Node i-esimo ya esta no se añade
      if (N_nodes.containsKey("N"+i)){}
      //caso contrario se añade a la lista
      else {N_nodes.put("N"+i,x);}
    }
    return N_nodes;
  }

//método Erdös Renyi para gener m aristas, devuelve objeto HashMap
//con parejas de nodos
//argumenmtos: Nodes(int), Edges(int), Directed(boolean), Auto(boolean)

  public static HashMap <String,Integer> genErdosRenyi(int n, int m, Boolean d, Boolean a){
    //creación de objeto HashMap para guardar la lista de aristas
    HashMap <String,Integer> G_ErdosRenyi = new HashMap <String,Integer>();
    //creación de objeto HashMap para guardar los Nodos
    HashMap <String,Double> Nodes = new HashMap <String,Double>();
    //variable para el número de nodos
    int edgs=m;
    int N_i;
    int N_f;
    //int num_nod=Nodes.size();
    //variable donde guardamos el objeto de la clase HashMap
    //creado por el metodo node que guarda los nodos
    Nodes = node(n);
    Boolean dir=d;
    Boolean Auto=a;
    for (int i=0;i<edgs;i++){
      double n_i=(Nodes.size())*Math.random();
      //System.out.println(n_i);
      double n_f=(Nodes.size())*Math.random();
      //System.out.println(n_f);
      N_i=(int)n_i;
      N_f=(int)n_f;
      //condició para evitar repeticiones en lista de aristas
      if (G_ErdosRenyi.containsKey("N"+N_i+"->N"+N_f)){}
      else{
        if (dir){
          if (Auto && (N_i == N_f)){ G_ErdosRenyi.put("N"+N_i+"->N"+N_f,i);System.out.println(N_i+"->"+N_f);}
          if (!Auto && N_i!=N_f){G_ErdosRenyi.put("N"+N_i+"->N"+N_f,i); System.out.println(N_i+"->"+N_f);}
        }
        if (!dir){
          if (Auto && (N_i == N_f)){ G_ErdosRenyi.put("N"+N_i+"->N"+N_f+"\t[arrowhead=none]",i); System.out.println(N_i+"--"+N_f);}
          if (!Auto && N_i!=N_f) {G_ErdosRenyi.put("N"+N_i+"->N"+N_f+"\t[arrowhead=none]",i); System.out.println(N_i+"--"+N_f);}
        }
      }
    }
    //regresa objeto de la clase HashMap con las aristas generadas por
    //genErdosRenyi
    return G_ErdosRenyi;
  }

//metodo Gilbert para generacion de m aristas para cada par del conjunto de
// nodos con probabilidad System.out.print();

  public static HashMap <String,Integer> genGilbert(int n, int m, double p, Boolean d, Boolean a){
    HashMap<String,Integer> G_Gilbert = new HashMap<String,Integer>();
    HashMap<String,Double> Nodos_glbr = new HashMap<String,Double>();
    Nodos_glbr = node(n);
    //variables locales
    int num_nod=Nodos_glbr.size();
    int edgs=m;
    double prbl=p;
    Boolean dir=d;
    Boolean Auto=a;
    int edgs_cnt=0;
    //ciclos para recorrer todas las parejas posibles de nodos y asignar arista
    //si numero aleatorio es superior a p
    for (int i=0; i<num_nod; i++){
      for(int j=0; j<num_nod; j++){
	      if (G_Gilbert.containsKey("N"+i+"->N"+j)){}
	      else{
		double rnd=Math.random();
		if (rnd>prbl){
		  if(dir){
		    if(Auto && i==j){G_Gilbert.put("N"+i+"->N"+j,edgs_cnt);}
		    if(!Auto && i!=j){G_Gilbert.put("N"+i+"->N"+j,edgs_cnt);}
		  }
		  if(!dir){
		    if(Auto && i==j){G_Gilbert.put("N"+i+"->N"+j+"\t[arrowhead=none]",edgs_cnt);}
		    if(!Auto && i!=j){G_Gilbert.put("N"+i+"->N"+j+"\t[arrowhead=none]",edgs_cnt);}
		  }
		}
		edgs_cnt=+edgs_cnt;
	}
      }
    }
    //regresa objeto de la clase HashMap con las aristas generadas por
    //genGilbert
    return G_Gilbert;
  }

  //metodo Geografico simple; regresa objeo tipo HashMap con las aristas generadas
  //por la regla de la r vecindad.
  public static HashMap <String,Integer> genSimpleGeo(int n, double r, Boolean d, Boolean a){
    HashMap<String,Integer> G_SimpleGeo= new HashMap<String,Integer>();
    HashMap <String,Double> N_x = new HashMap <String,Double>();
    HashMap <String,Double> N_y = new HashMap <String,Double>();
 
 

    int N=n;
    int edgs_cnt;
    double B=r;
    double R=0;
    Boolean Auto=a;
    Boolean dir=d;
    N_x=node(N);
    N_y=node(N);

    double x1=0;
    double x2=0;
    double y1=0;
    double y2=0;    
    

    edgs_cnt=0;

    for (int i=0; i<N; i++){
       for (int j=0; j<N; j++){
	  x1=N_x.get("N"+i);
	  x2=N_x.get("N"+j);
	  x1=N_y.get("N"+i);
	  y2=N_y.get("N"+j);
	  R=Math.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));	
          if(R<B){
	      if (G_SimpleGeo.containsKey("N"+i+"->N"+j)){}
	      else{
		  if(dir){
		    if(Auto && i==j){G_SimpleGeo.put("N"+i+"->N"+j,edgs_cnt);}
		    if(!Auto && i!=j){G_SimpleGeo.put("N"+i+"->N"+j,edgs_cnt);}
		  }
		  if(!dir){
		    if(Auto && i==j){G_SimpleGeo.put("N"+i+"->N"+j+"\t [arrowhead=none]",edgs_cnt);}
		    if(!Auto && i!=j){G_SimpleGeo.put("N"+i+"->N"+j+"\t [arrowhead=none]",edgs_cnt);}
		  }
		}
	  }
	 edgs_cnt=edgs_cnt+1;
       }
    }

    return G_SimpleGeo;
  }

  //metodo Barabasi Albert; regresa objeto tipo HashMap con las aristas generadas con
  //la probabilidad p=1-deg(nodo)/d.
  public static HashMap <String,Integer> genBarabasiAlbert(int n, int D_p, Boolean d, Boolean a){

	HashMap<String,Integer> G_BarabasiAlbert= new HashMap<String,Integer>();
	HashMap<String,Double> nodos = new HashMap<String,Double>();

	int N=n;
	int D=D_p;
	int edgs_cnt=0;
	double p=0;
	double deg_v=0;
	double rand=0;
	Boolean dir=d;
	Boolean Auto=a;
	nodos=node(N);
	
	
	for (int h=0; h<N;h++){
		nodos.put("N"+h,0.0);
	}	

	for (int i=0; i<N; i++){	
		for (int j=i; j<N; j++){
		      if (G_BarabasiAlbert.containsKey("N"+i+"->N"+j)){}
      		      else{
			
				deg_v=nodos.get("N"+i);
				p=1.0-deg_v/(double)D;
				rand=Math.random();
				if(rand < p){
					if(dir){
						if (Auto && ("N"+i=="N"+j)){
							G_BarabasiAlbert.put("N"+i+"->N"+j,edgs_cnt);
							double temp_deg=nodos.get("N"+i);				
							nodos.put("N"+i,temp_deg+1);
							temp_deg=nodos.get("N"+j);
							nodos.put("N"+j,temp_deg+1);
							edgs_cnt=edgs_cnt+1;
						}
						if (!Auto && ("N"+i!="N"+j)){
							G_BarabasiAlbert.put("N"+i+"->N"+j,edgs_cnt);
							double temp_deg=nodos.get("N"+i);				
							nodos.put("N"+i,temp_deg+1);
							temp_deg=nodos.get("N"+j);
							nodos.put("N"+j,temp_deg+1);
							edgs_cnt=edgs_cnt+1;
						}
					}
					if(!dir){
						if (Auto && ("N"+i=="N"+j)){
							G_BarabasiAlbert.put("N"+i+"->N"+j+"\t [arrowhead=none]",edgs_cnt);
							double temp_deg=nodos.get("N"+i);				
							nodos.put("N"+i,temp_deg+1);
							temp_deg=nodos.get("N"+j);
							nodos.put("N"+j,temp_deg+1);
							edgs_cnt=edgs_cnt+1;
						}
						if (!Auto && ("N"+i!="N"+j)){
							G_BarabasiAlbert.put("N"+i+"->N"+j+"\t [arrowhead=none]",edgs_cnt);
							double temp_deg=nodos.get("N"+i);				
							nodos.put("N"+i,temp_deg+1);
							temp_deg=nodos.get("N"+j);
							nodos.put("N"+j,temp_deg+1);
							edgs_cnt=edgs_cnt+1;
						}
					}
				}
			}
		}
		
	}	

    return G_BarabasiAlbert;
  }

  //metodo para guardar la lista generada por los diferentes metodos en archivo .viz
  public void toViz(String file_name, Boolean d, HashMap<String,Integer> Graph_method) throws IOException{

      Boolean dir=d;
      File file =new File(file_name+".viz");
      file.createNewFile();
      FileWriter writer = new FileWriter(file);
      
      writer.write("Digraph G {\n");
      Graph_method.forEach((k,v) -> {
            try {
		writer.write("\t"+k+	";\n");
	    } catch(IOException e) {
                System.out.println("An error occurred.");
                e.printStackTrace();
            }
            
       });

      writer.write("}");
      writer.flush();
      writer.close();
  }


  public static void main(String[] args){
    // creacion de objetos de la clase Graph
    Graph g1 = new Graph();
    Graph g2 = new Graph();
    Graph g3 = new Graph();
    Graph g4 = new Graph();

    //variable numero de nodos
    int n=0;
    //variable numero de aristas
    int m=0;
    //variable metodo Gilbert
    double p=0.9973;
    //variable metodo Geografica Simple
    double r=0.08;
    //variable metodo Barabasi Albert
    int D=5;
    //variable Dirigido
    Boolean d=false;
    //variable Autoconectado
    Boolean a=false;
    //objetos de la clase HashMap para guardar el regreso de los metodos
    HashMap <String,Double> Nodes = new HashMap <String,Double>();
    HashMap <String,Integer> ErdRny = new HashMap <String,Integer>();
    HashMap <String,Integer> Gilbert = new HashMap <String,Integer>();
    HashMap <String,Integer> SimpleGeo = new HashMap <String,Integer>();
    HashMap <String,Integer> BarabasiAlbert = new HashMap <String,Integer>();

    Scanner keyboard = new Scanner(System.in);


    System.out.println("Numero de nodos");
    n=keyboard.nextInt();

    System.out.println("Numero de aristas ");
    m=keyboard.nextInt();
/*
    System.out.println("Probabilidad para metodo Gilbert");
    p=keyboard.nextDouble();

    System.out.println("Radio de la vecindad para metodo Geografico simple");
    r=keyboard.nextDouble();

    System.out.println("Parametro d para metodo Barabasi Albert");
    D=keyboard.nextInt();

    System.out.println("Dirigido");
    d=keyboard.nextBoolean();

    System.out.println("Autoconectado");
    a=keyboard.nextBoolean();
*/
			
    Nodes= g1.node(n);

    ErdRny=g1.genErdosRenyi(n,m,d,a);
    Gilbert=g2.genGilbert(n,m,p,d,a);
    SimpleGeo=g3.genSimpleGeo(n,r,d,a);
    BarabasiAlbert=g4.genBarabasiAlbert(n,D,d,a);

	try {
	    g1.toViz("ErdRny_n"+n+"_m"+m, d, ErdRny);
	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
	}

	try {
	    g1.toViz("Gilbert_n"+n+"_p"+p, d, Gilbert);
	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
	}

	try {
	    g1.toViz("SimpleGeo_n"+n+"_r"+r, d, SimpleGeo);
	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
	}

	try {
	    g1.toViz("BarabasiAlbert_n"+n+"_d"+D, d, BarabasiAlbert);
	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
	}

    keyboard.close();
  }

}
