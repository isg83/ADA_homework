import java.util.*;
import java.io.*;
import java.nio.file.DirectoryIteratorException;

public class Grafo{
    HashMap<String, Double> Nodes = new HashMap<String,Double>();
    HashMap<String, Double> Edges = new HashMap<String,Double>();
    
    public Grafo(HashMap<String, Double> Nodes, HashMap<String, Double> Edges){   
        this.Nodes=Nodes;
        this.Edges=Edges;
    }

    public void nodes(Integer N){
        int n=N;
        for (int i=0;i<=n;i++){
            this.Nodes.put("N"+i,0.0);            
        }
    }

    public void ER(Integer m, Boolean Dir, Boolean Auto){
    
        Integer M=m;
        Boolean D=Dir;
        Boolean A=Auto;
        Double p;
    
        for (int j=0; j<=M; j++){
            int Ni=(int)(this.Nodes.size()*Math.random()+1.0);
            int Nf=(int)(this.Nodes.size()*Math.random()+1.0);

            if (D){
                if (this.Edges.containsKey("N"+Ni+"->"+"N"+Nf)){}
                else{this.Edges.put("N"+Ni+"->"+"N"+Nf,1.0);}
            }
            if(!D){
                if (this.Edges.containsKey("N"+Ni+"->"+"N"+Nf) || this.Edges.containsKey("N"+Nf+"->"+"N"+Ni)){}
                else{this.Edges.put("N"+Ni+"--"+"N"+Nf,1.0);}
            }
        }
    }
    public static void main(String[] args){
        HashMap<String,Double> N= new HashMap<>();
        HashMap<String,Double> E= new HashMap<>();

        int n=400;        
        int m=400;
        Boolean d=false;
        Boolean a=false;
        
        Grafo G1= new Grafo(N,E);
        G1.nodes(n);
        G1.ER(m,d,a);

        for(String key: G1.Edges.keySet()){
            System.out.println(key);
        }
    }
}