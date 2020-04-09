import java.util.*;
import java.io.*;
import java.nio.file.DirectoryIteratorException;

public class Small_DFS_recursive{

    public HashMap<String, Integer> DFS_HM = new HashMap<>();
    public Map<String, Integer> Edges;
    public Map<String,Double> Nodos_Graph;    
    public Map<String,Boolean> NotDiscover_map;
    public Stack<String> stack;
    public String Nodo_previo;  

    public Small_DFS_recursive(Map<String, Integer> Edges,HashMap<String, Integer> DFS_HM, Map<String,Double> Nodos_Graph,String Nodo_previo, Map<String,Boolean> NotDiscover_map,Stack<String> stack){
        System.out.println("creando objeto de clase small");
        this.Edges=Edges;
        this.DFS_HM=DFS_HM;
        this.Nodos_Graph=Nodos_Graph;
        this.Nodo_previo=Nodo_previo;
        this.NotDiscover_map=NotDiscover_map;
        this.stack=stack;
        /*
        for (String key0: this.Nodos_Graph.keySet()){
            System.out.println(key0);
        }
        */
        //System.out.println(this.Nodo_previo);
    }
//pero regreso al mismo problema en esta función debo decirle que tipo de retorno da el método
//lo que quiero que regrese es el objeto modificado. Si le pongo void no se como terminar 
//la recursividad pues no puedo return null  
    public  void small_DFSr(String N_prev){
      
        String key1=N_prev;      
        this.NotDiscover_map.put(key1,false);      
        //System.out.println(key1+","+this.NotDiscover_map.get(key1));
        //System.out.println(this.Nodos_Graph.size());
        //System.out.println(this.stack.size());
        for(String key2: this.Nodos_Graph.keySet()){   
          //System.out.println(this.Edges.containsKey(key1+"--"+key2)+","+this.NotDiscover_map.get(key2));           
          if(this.Edges.containsKey(key1+"--"+key2) && this.NotDiscover_map.get(key2)){            
            this.DFS_HM.put(key1+"--"+key2,1);
            //System.out.println(key1+"--"+key2);
            this.NotDiscover_map.put(key2,false);
            this.stack.push(key2);
          }       
        }     
        
        try{
            key1=this.stack.pop();
            small_DFSr(key1);
        }
        catch (EmptyStackException e) {
            System.out.println("se alcanzó el final del stack");
        }             
        
     
    }
    
    public static void main(String[] args){
        //Small_DFS_recursive object= new Small_DFS_recursive(Graph, DFS_HM, Nodos_Graph, Nodo_previo, NotDiscover_map, stack);
        //object.small_DFSr(N_prev);
    }
    
}