import java.util.*;
import java.io.*;
import java.nio.file.DirectoryIteratorException;

public class Small_DFS_recursive{

    public HashMap<String, Integer> Graph;
    public HashMap<String, Integer> DFS_HM;
    public Map<String,Integer> Nodos_Graph;
    public String Nodo_previo;
    public Map<String,Boolean> NotDiscover_map;
    public Stack<String> stack;

    public Small_DFS_recursive(HashMap<String, Integer> Graph,HashMap<String, Integer> DFS_HM, Map<String,Integer> Nodos_Graph,String Nodo_previo, Map<String,Boolean> NotDiscover_map,Stack<String> stack){

        this.Graph=Graph;
        this.DFS_HM=DFS_HM;
        this.Nodos_Graph=Nodos_Graph;
        this.Nodo_previo=Nodo_previo;
        this.NotDiscover_map=NotDiscover_map;
        this.stack=stack;
    }
//pero regreso al mismo problema en esta función debo decirle que tipo de retorno da el método
//lo que quiero que regrese es el objeto modificado. Si le pongo void no se como terminar 
//la recursividad pues no puedo return null  
    public  void small_DFSr(String N_prev){
      
        String key1=N_prev;      
        this.NotDiscover_map.put(key1,false);      
        for(String key2: this.Nodos_Graph.keySet()){
          if(this.Graph.containsKey(key1+"--"+key2) && this.NotDiscover_map.get(key2)){
            this.DFS_HM.put(key1+"--"+key2,0);
            this.NotDiscover_map.put(key2,false);
            this.stack.push(key2);
          }            
        }
        try{
            key1=this.stack.pop();
            small_DFSr(key1);
        }
        catch (EmptyStackException e) {
            break;
        }  
        
    }
    
    public static void main(String[] args){
        Small_DFS_recursive object= new Small_DFS_recursive(Graph, DFS_HM, Nodos_Graph, Nodo_previo, NotDiscover_map, stack);
        object.small_DFSr(N_prev);

    }
    
}