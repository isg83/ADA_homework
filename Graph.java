import java.util.*;
import java.util.stream.Collectors;
//import sun.misc.Queue;
import java.io.*;
import java.nio.file.DirectoryIteratorException;
import java.text.DecimalFormat;

public class Graph{
//metodo devuelve objeto HashMap con N nodos, argumento N
  public static HashMap <String,Double> node(int N) {
    // creacion de objeto de la clase HashMap llamado N_nodes
    HashMap<String, Double> N_nodes = new HashMap<String, Double>();
    // variable de metodo para contener el número de Nodos
    int nodes = N;
    double x = 0;
    // ciclo para meter nodo a la lista HashMap
    for (int i = 0; i < nodes; i++) {
      x = Math.random();
      // si el Node i-esimo ya esta no se añade
      if (N_nodes.containsKey("N" + i)) {
      }
      // caso contrario se añade a la lista
      else {
        N_nodes.put("N" + i, x);
      }
    }
    return N_nodes;
  }

  // método Erdös Renyi para gener m aristas, devuelve objeto HashMap
  // con parejas de nodos
  // argumenmtos: Nodes(int), Edges(int), Directed(boolean), Auto(boolean)

  public static HashMap<String, Integer> genErdosRenyi(int n, int m, Boolean d, Boolean a) {
    // creación de objeto HashMap para guardar la lista de aristas
    HashMap<String, Integer> G_ErdosRenyi = new HashMap<String, Integer>();
    // creación de objeto HashMap para guardar los Nodos
    HashMap<String, Double> Nodes = new HashMap<String, Double>();
    // variable para el número de nodos
    int edgs = m;
    int N_i;
    int N_f;
    int w;
    // int num_nod=Nodes.size();
    // variable donde guardamos el objeto de la clase HashMap
    // creado por el metodo node que guarda los nodos
    Nodes = node(n);
    Boolean dir = d;
    Boolean Auto = a;
    for (int i = 0; i < edgs; i++) {
      double n_i = (Nodes.size()) * Math.random();
      // System.out.println(n_i);
      double n_f = (Nodes.size()) * Math.random();
      // System.out.println(n_f);
      N_i = (int) n_i;
      N_f = (int) n_f;
      // condición para evitar repeticiones en lista de aristas
      if (G_ErdosRenyi.containsKey("N" + N_i + "->N" + N_f)) { } 
      else {
        double w_d=50*Math.random()+1;
        w=(int) w_d;
        if (dir) {
          if (Auto && (N_i == N_f)) {            
            G_ErdosRenyi.put("N" + N_i + "->N" + N_f, w);
            //System.out.println(N_i + "->" + N_f);
          }
          if (!Auto && N_i != N_f) {
            G_ErdosRenyi.put("N" + N_i + "->N" + N_f, w);
            //System.out.println(N_i + "->" + N_f);
          }
        }
        if (!dir && !G_ErdosRenyi.containsKey("N" + N_f + "--N" + N_i)) {
          if (Auto && (N_i == N_f)) {
            G_ErdosRenyi.put("N" + N_i + "--N" + N_f , w);
            //System.out.println(N_i + "--" + N_f);
          }
          if (!Auto && N_i != N_f) {
            G_ErdosRenyi.put("N" + N_i + "--N" + N_f , w);
            //System.out.println(N_i + "--" + N_f);
          }
        }
      }
    }
    // regresa objeto de la clase HashMap con las aristas generadas por
    // genErdosRenyi
    return G_ErdosRenyi;
  }

  // metodo Gilbert para generacion de m aristas para cada par del conjunto de
  // nodos con probabilidad System.out.print();

  public static HashMap<String, Integer> genGilbert(int n, int m, double p, Boolean d, Boolean a) {
    HashMap<String, Integer> G_Gilbert = new HashMap<String, Integer>();
    HashMap<String, Double> Nodos_glbr = new HashMap<String, Double>();
    Nodos_glbr = node(n);
    // variables locales
    int num_nod = Nodos_glbr.size();
    int edgs = m;
    double prbl = p;
    Boolean dir = d;
    Boolean Auto = a;
    //int edgs_wght = 1;
    // ciclos para recorrer todas las parejas posibles de nodos y asignar arista
    // si numero aleatorio es superior a p
    for (int i = 0; i < num_nod; i++) {
      for (int j = 0; j < num_nod; j++) {
        if (G_Gilbert.containsKey("N" + i + "->N" + j)) { } 
        else {
          double rnd = Math.random();
          double edgs_wght_d=100*Math.random()+1.0;
          int edgs_wght=(int) edgs_wght_d;
          if (rnd > prbl) {
            if (dir) {
              if (Auto && i == j) {
                G_Gilbert.put("N" + i + "->N" + j, edgs_wght);
              }
              if (!Auto && i != j) {
                G_Gilbert.put("N" + i + "->N" + j, edgs_wght);
              }
            }
            if (!dir && !G_Gilbert.containsKey("N" + j + "--N" + i )) {
              if (Auto && i == j) {
                G_Gilbert.put("N" + i + "--N" + j , edgs_wght);
              }
              if (!Auto && i != j) {
                G_Gilbert.put("N" + i + "--N" + j , edgs_wght);
              }
            }
          }
        }
      }
    }
    // regresa objeto de la clase HashMap con las aristas generadas por
    // genGilbert
    return G_Gilbert;
  }
//arreglar la duplicidad en no dirigidos
  // metodo Geografico simple; regresa objeo tipo HashMap con las aristas
  // generadas
  // por la regla de la r vecindad.
  public static HashMap<String, Integer> genSimpleGeo(int n, double r, Boolean d, Boolean a) {

    HashMap<String, Integer> G_SimpleGeo = new HashMap<String, Integer>();
    HashMap<String, Double> N_x = new HashMap<String, Double>();
    HashMap<String, Double> N_y = new HashMap<String, Double>();

    int N = n;
    int edgs_cnt;
    double B = r;
    double R = 0;
    Boolean Auto = a;
    Boolean dir = d;
    N_x = node(N);
    N_y = node(N);

    double x1 = 0;
    double x2 = 0;
    double y1 = 0;
    double y2 = 0;

    edgs_cnt = 0;

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        x1 = N_x.get("N" + i);
        x2 = N_x.get("N" + j);
        x1 = N_y.get("N" + i);
        y2 = N_y.get("N" + j);
        R = Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        if (R < B) {
          if (G_SimpleGeo.containsKey("N" + i + "->N" + j)) { } 
          else {
            double edgs_wght_d=100*Math.random()+1.0;
            int edgs_wght=(int) edgs_wght_d;
            if (dir) {
              if (Auto && i == j) {
                G_SimpleGeo.put("N" + i + "->N" + j, edgs_wght);
              }
              if (!Auto && i != j) {
                G_SimpleGeo.put("N" + i + "->N" + j, edgs_wght);
              }
            }
            if (!dir && !G_SimpleGeo.containsKey("N" + j + "--N" + i)) {
              if (Auto && i == j) {
                G_SimpleGeo.put("N" + i + "--N" + j , edgs_wght);
              }
              if (!Auto && i != j) {
                G_SimpleGeo.put("N" + i + "--N" + j , edgs_wght);
              }
            }
          }
        }
        //edgs_cnt = edgs_cnt + 1;
      }
    }

    return G_SimpleGeo;
  }

  // metodo Barabasi Albert; regresa objeto tipo HashMap con las aristas generadas
  // con
  // la probabilidad p=1-deg(nodo)/d.
  public static HashMap<String, Integer> genBarabasiAlbert(final int n, final int D_p, final Boolean d,
      final Boolean a) {

    final HashMap<String, Integer> G_BarabasiAlbert = new HashMap<String, Integer>();
    HashMap<String, Double> nodos = new HashMap<String, Double>();

    final int N = n;
    final int D = D_p;
    int edgs_cnt = 0;
    double p = 0;
    double deg_v = 0;
    double rand = 0;
    final Boolean dir = d;
    final Boolean Auto = a;
    nodos = node(N);

    for (int h = 0; h < N; h++) {
      nodos.put("N" + h, 0.0);
    }

    for (int i = 0; i < N; i++) {
      for (int j = i; j < N; j++) {
        if (G_BarabasiAlbert.containsKey("N" + i + "->N" + j)) {
        } else {

          deg_v = nodos.get("N" + i);
          p = 1.0 - deg_v / (double) D;
          rand = Math.random();
          if (rand < p) {
            double edgs_wght_d=100*Math.random()+1.0;
            int edgs_wght=(int) edgs_wght_d;
            if (dir) {
              if (Auto && ("N" + i == "N" + j)) {
                G_BarabasiAlbert.put("N" + i + "->N" + j, edgs_wght);
                double temp_deg = nodos.get("N" + i);
                nodos.put("N" + i, temp_deg + 1);
                temp_deg = nodos.get("N" + j);
                nodos.put("N" + j, temp_deg + 1);
                //edgs_cnt = edgs_cnt + 1;
              }
              if (!Auto && ("N" + i != "N" + j)) {
                G_BarabasiAlbert.put("N" + i + "->N" + j, edgs_wght);
                double temp_deg = nodos.get("N" + i);
                nodos.put("N" + i, temp_deg + 1);
                temp_deg = nodos.get("N" + j);
                nodos.put("N" + j, temp_deg + 1);
                //edgs_cnt = edgs_cnt + 1;
              }
            }
            if (!dir && !G_BarabasiAlbert.containsKey("N" + j + "--N" + i)) {
              if (Auto && ("N" + i == "N" + j)) {
                G_BarabasiAlbert.put("N" + i + "--N" + j , edgs_wght);
                double temp_deg = nodos.get("N" + i);
                nodos.put("N" + i, temp_deg + 1);
                temp_deg = nodos.get("N" + j);
                nodos.put("N" + j, temp_deg + 1);
                //edgs_cnt = edgs_cnt + 1;
              }
              if (!Auto && ("N" + i != "N" + j)) {
                G_BarabasiAlbert.put("N" + i + "--N" + j , edgs_wght);
                double temp_deg = nodos.get("N" + i);    
                nodos.put("N" + i, temp_deg + 1);
                temp_deg = nodos.get("N" + j);
                nodos.put("N" + j, temp_deg + 1);
                //edgs_cnt = edgs_cnt + 1;
              }
            }
          }
        }
      }

    }

    return G_BarabasiAlbert;
  }

  // metodo para guardar la lista generada por los diferentes metodos en archivo
  // .dot en un principio generaba archivos viz pero gephi no los puede leer.
  //mantengo el nombre pues es el indicado en las instrucciones dadas en clase.
  public void toViz(final String file_name, final Boolean d, final HashMap<String, Integer> Graph_method)
      throws IOException {

    final Boolean dir = d;
    final File file = new File(file_name + ".dot");
    file.createNewFile();
    final FileWriter writer = new FileWriter(file);
    
    if(dir){
    writer.write("Digraph G {\n");
    Graph_method.forEach((k, v) -> {
      try {
        writer.write("\t" + k + "[weight=" + v + "];\n");
      } catch (final IOException e) {
        System.out.println("An error occurred.");
        e.printStackTrace();
      }
    });
    }

    if(!dir){
      writer.write("Graph G {\n");
      Graph_method.forEach((k, v) -> {
        try {
          writer.write("\t" + k + "[weight=" + v + "];\n");
        } catch (final IOException e) {
          System.out.println("An error occurred.");
          e.printStackTrace();
        }
      });
      }
      

    writer.write("}");
    writer.flush();
    writer.close();
  }

   // metodo para guardar la lista generada por los diferentes metodos en archivo
  // .dot para el caso donde el Map contiene double 
  public void toVizDoub(String file_name, Boolean d, Map<String, Double> Graph_method)
      throws IOException {

    Boolean dir = d;
    File file = new File(file_name + ".dot");
    file.createNewFile();
    FileWriter writer = new FileWriter(file);
    
    if(dir){
    writer.write("Digraph G {\n");
    Graph_method.forEach((k, v) -> {
      try {
        writer.write("\t" + k + "[weight=" + v + "];\n");
      } catch (final IOException e) {
        System.out.println("An error occurred.");
        e.printStackTrace();
      }
    });
    }
DecimalFormat df = new DecimalFormat("#.00");

    if(!dir){
      writer.write("Graph G {\n");
      Graph_method.forEach((k, v) -> {
        try {
          
          writer.write("\t" + k + "\t [label="+ k +"(" + df.format(v) + ")];\n");
        } catch (final IOException e) {
          System.out.println("An error occurred.");
          e.printStackTrace();
        }
      });
    }
    writer.write("}");
    writer.flush();
    writer.close();
  }

//metodo BFS de forma iterativa.
//en un principio no encontraba como iterar un Hashmap pues los ejemplos que buscaba siempre usaban a map 
//en lugar de un Hashmap pero ya se que es posible. 

 public HashMap<String,Integer> graphToBFS(final HashMap<String, Integer> Graph_in,HashMap<String,Double> Nodos_G_in, int N_source,Boolean Dir){
  HashMap<String,Integer> G_BFS = new HashMap<String,Integer>();
  Map<String,Integer> map_to_hold_the_bloody_and_uniterable_HashMap_edges = new HashMap<String,Integer>(Graph_in);
  Map<String,Double> another_map_to_hold_the_other_HashMap_nodes = new HashMap<String,Double>(Nodos_G_in);
  Map<String,Integer> layer_map = new HashMap<String,Integer>();
  Map<String,Boolean> NotDiscover_map = new HashMap<String,Boolean>();

  int Total_nodos_grafo=Nodos_G_in.size();
  int Nodo_previo=N_source;
  int edges_counter=0;
  int l=0;
  int Numero_Nodos_descubiertos=1;
  int nodos_restantes;
  int detached=0;
  Boolean Dirigido=Dir;

  for (Map.Entry<String, Double> entry : another_map_to_hold_the_other_HashMap_nodes.entrySet()) {
    String key = entry.getKey();
    Object value = entry.getValue();
    layer_map.put(key,Total_nodos_grafo);
    NotDiscover_map.put(key, true);
  }
  
  layer_map.put("N"+Nodo_previo, 0);
  NotDiscover_map.put("N"+Nodo_previo, false);

if(!Dirigido){  
  while(Numero_Nodos_descubiertos<Total_nodos_grafo-1){ 
    nodos_restantes=Total_nodos_grafo-Numero_Nodos_descubiertos;  
    //System.out.println(Numero_Nodos_descubiertos+","+Total_nodos_grafo); 
    for(String key1: layer_map.keySet()){      
      if(layer_map.get(key1)==l){
        //System.out.println(layer_map.get(key1));
        for(String key2: another_map_to_hold_the_other_HashMap_nodes.keySet()){
          
          if(map_to_hold_the_bloody_and_uniterable_HashMap_edges.containsKey(key1+"--"+key2) && NotDiscover_map.get(key2)){
            //System.out.println(key1+","+key2);
            G_BFS.put(key1+"--"+key2, map_to_hold_the_bloody_and_uniterable_HashMap_edges.get(key1+"--"+key2));
            NotDiscover_map.put(key2, false);
            layer_map.put(key2,l+1);
            Numero_Nodos_descubiertos=Numero_Nodos_descubiertos+1;
            edges_counter=edges_counter+1;
          }
        }
      }
    }
    l=l+1;
    int temp=nodos_restantes;
    nodos_restantes=Total_nodos_grafo-Numero_Nodos_descubiertos;
    if (nodos_restantes==temp){
      detached=detached+1;
      if (detached>2){
        //System.out.println("Desconectada");
        break;
      }
    }
    
  }
}

if(Dirigido){
  while(Numero_Nodos_descubiertos<Total_nodos_grafo-1){
    nodos_restantes=Total_nodos_grafo-Numero_Nodos_descubiertos;  
    for(String key1: layer_map.keySet()){
      if(layer_map.get(key1)==l){
        for(String key2: another_map_to_hold_the_other_HashMap_nodes.keySet()){
          if(map_to_hold_the_bloody_and_uniterable_HashMap_edges.containsKey(key1+"->"+key2) && NotDiscover_map.get(key2)){
            G_BFS.put(key1+"->"+key2, map_to_hold_the_bloody_and_uniterable_HashMap_edges.get(key1+"->"+key2));
            NotDiscover_map.put(key2, false);
            layer_map.put(key2,l+1);
            Numero_Nodos_descubiertos=Numero_Nodos_descubiertos+1;
            edges_counter=edges_counter+1;
          }
        }
      }
    }
    l=l+1;
    int temp=nodos_restantes;
    nodos_restantes=Total_nodos_grafo-Numero_Nodos_descubiertos;
    if (nodos_restantes==temp){
      detached=detached+1;
      if (detached>2){break;}
    }
  }
}

  return G_BFS;
 }

//Método para obtener el DFS de forma iterativa.

public HashMap<String, Integer> graphToDFSi(HashMap<String, Integer> Graph_in,HashMap<String,Double> Nodos_G_in, int N_source,Boolean Dir){
  HashMap<String,Integer> G_DFSi = new HashMap<String,Integer>();
  Map<String,Integer> map_to_hold_the_bloody_and_uniterable_HashMap_edges = new HashMap<String,Integer>(Graph_in);
  Map<String,Double> another_map_to_hold_the_other_HashMap_nodes = new HashMap<String,Double>(Nodos_G_in);
  Stack<String> stack = new Stack<String>();
  Map<String,Boolean> NotDiscover_map = new HashMap<String,Boolean>();

  int Total_nodos_grafo=Nodos_G_in.size();
  int Nodo_previo=N_source;
  int Numero_Nodos_descubiertos=1;
  int nodos_restantes;
  String key2;
  Boolean Dirigido=Dir;

  for(String key1: another_map_to_hold_the_other_HashMap_nodes.keySet()){
    NotDiscover_map.put(key1, true);
  }
    
  NotDiscover_map.put("N"+Nodo_previo, false);

  stack.push("N"+Nodo_previo);
  key2=stack.pop();

  if (Dir){
    while(true){    
      for (String key1: another_map_to_hold_the_other_HashMap_nodes.keySet()){
        if(map_to_hold_the_bloody_and_uniterable_HashMap_edges.containsKey(key2+"->"+key1) && NotDiscover_map.get(key1)){
          G_DFSi.put(key2+"->"+key1, map_to_hold_the_bloody_and_uniterable_HashMap_edges.get(key2+"->"+key1));
          stack.push(key1);  
          NotDiscover_map.put(key1, false);          
        }      
      }
      try{
        key2=stack.pop();
      }
      catch (EmptyStackException e) {
        break;
     }  
    }
  }

  if (!Dir){
    while(true){    
      for (String key1: another_map_to_hold_the_other_HashMap_nodes.keySet()){
        if(map_to_hold_the_bloody_and_uniterable_HashMap_edges.containsKey(key2+"--"+key1) && NotDiscover_map.get(key1)){
          G_DFSi.put(key2+"--"+key1, map_to_hold_the_bloody_and_uniterable_HashMap_edges.get(key2+"--"+key1));
          stack.push(key1);  
          NotDiscover_map.put(key1, false);          
        }      
      }
      try{
        key2=stack.pop();
      }
      catch (EmptyStackException e) {
        break;
     }  
    }
  }

  return G_DFSi;
}

//método para obtener DFS de forma recursiva

public HashMap<String, Integer> graphToDFSr(HashMap<String, Integer> Graph_in,HashMap<String,Double> Nodos_G_in, int N_source,Boolean Dir){
  HashMap<String,Integer> G_DFSr = new HashMap<String,Integer>();
  HashMap<String,Integer> G_DFSrO = new HashMap<String,Integer>();
  Map<String,Integer> map_to_hold_the_bloody_and_uniterable_HashMap_edges = new HashMap<String,Integer>(Graph_in);
  Map<String,Double> another_map_to_hold_the_other_HashMap_nodes = new HashMap<String,Double>(Nodos_G_in);
  Map<String,Boolean> NotDiscover_map = new HashMap<String,Boolean>();
  Stack<String> stack = new Stack<String>();

  String Nodo_previo="N"+N_source;
  String key1;
  Boolean Dirigido=Dir;

  for (String key: another_map_to_hold_the_other_HashMap_nodes.keySet()){
    NotDiscover_map.put(key,true);
    //System.out.println(key+"true");
  }

  Small_DFS_recursive hero = new Small_DFS_recursive(map_to_hold_the_bloody_and_uniterable_HashMap_edges, G_DFSr, another_map_to_hold_the_other_HashMap_nodes, Nodo_previo, NotDiscover_map, stack);

    hero.small_DFSr(Nodo_previo);  

  //System.out.println("imprimiendo salida de recursividad en atributo DFS");
  /*
  for(String keyO: hero. DFS_HM.keySet()){
    System.out.println(keyO);
  }
  */
  G_DFSrO.putAll(hero.DFS_HM);
  

  return G_DFSrO;
}

//metodo para ordenar map por valor

public static Map<String,Double> sortByValue(final Map<String,Double> Dist_al_nodo){
  return Dist_al_nodo.entrySet()
          .stream()
          .sorted((Map.Entry.<String, Double>comparingByValue()))
          .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1,e2) -> e1, LinkedHashMap::new));
}
//metodo Dijkstra

public Map<String,Double> getDijkstra(HashMap<String, Double> Graph_in,HashMap<String,Double> Nodos_G_in, String N_source,Boolean Dir){
  Map<String,Double> Graph = new HashMap<String,Double>(Graph_in);
  Map<String,Double> map_Distancia_Nodos_G = new HashMap<String,Double>(Nodos_G_in);
  Map<String,Double> order_by_dist_nodos_G = new HashMap<String,Double>();
  Map<String,Boolean> InSetQ = new HashMap<String,Boolean>();
  Map<String,String> prev = new HashMap<String,String>();
  Map<String,Double> Graph_min =new HashMap<String,Double>();
  List<String> removedList = new ArrayList<String>();

  double inf=Double.POSITIVE_INFINITY;
  double dist;
  int not_new_counter;
  int Total_nodos=map_Distancia_Nodos_G.size();
  int j=0;
  String key0=N_source;

  for(String key1 : map_Distancia_Nodos_G.keySet()){
    map_Distancia_Nodos_G.put(key1,inf);    
    InSetQ.put(key1, true);
  }
  map_Distancia_Nodos_G.put(key0,0.0);
  
  not_new_counter=map_Distancia_Nodos_G.size();
  while(InSetQ.containsValue(true)){
    removedList.add(key0);
    for(String key1 : map_Distancia_Nodos_G.keySet()){
      if(Graph.containsKey(key0+"--"+key1) && InSetQ.get( key0) && InSetQ.get(key1) ){
        //System.out.println(key0+"--"+key1+","+InSetQ.get(key0));
        dist=Graph.get(key0+"--"+key1)+map_Distancia_Nodos_G.get(key0);        
                        
        if(dist < map_Distancia_Nodos_G.get(key1)){
          map_Distancia_Nodos_G.put(key1, dist);
          //Graph_min.put(key0+"--"+key1,dist);
          prev.put(key1, key0);
        }
      }
    }
    InSetQ.put(key0,false);
    
    //se re obtiene los nodos totales ordenados de menor a mayor peso
    order_by_dist_nodos_G=sortByValue(map_Distancia_Nodos_G);
    /*    
    for (String key2: order_by_dist_nodos_G.keySet()){
      System.out.println(key2+","+order_by_dist_nodos_G.get(key2));
    }
    */
        
    for (int i=0; i < removedList.size(); i++){
      j=i;
      order_by_dist_nodos_G.remove(removedList.get(i));        
    }
   
    try {
      Map.Entry<String,Double> entry = order_by_dist_nodos_G.entrySet().stream().findFirst().get();
      key0=entry.getKey();
    } catch (Exception e) {
      System.out.println("trataste de obtener elemento de map vacio");
      break;
    }

    int temp= order_by_dist_nodos_G.size();    
    if (not_new_counter==temp){
      System.out.println("Grafo desconectado; al calcular Dijkstra");
      break;
    }
    not_new_counter=temp;
  } 
  for(String key3: prev.keySet()){    
    Graph_min.put(prev.get(key3)+"--"+key3,map_Distancia_Nodos_G.get(key3));
  }

  for(String key4: map_Distancia_Nodos_G.keySet()){
    Graph_min.put(key4,map_Distancia_Nodos_G.get(key4));
  } 
return Graph_min;
}


public HashMap<String, Double> RandomEdgeValues(HashMap<String, Integer> Graph_in, Double min, Double max){
  
  double mn=min;
  double mx=max;
  HashMap<String, Integer> Graph_to_hold = new HashMap<String, Integer>(Graph_in);
  HashMap<String, Double> Graph_to_return = new HashMap<String, Double>();

  for(String key0: Graph_to_hold.keySet()){
    double rnd_wgt=mn+(mx-mn)*Math.random(); 
    Graph_to_return.put(key0,rnd_wgt);
    //System.out.println(key0+" peso "+rnd_wgt);
  }

  return Graph_to_return;
}

public static void main(final String[] args){
    // creacion de objetos de la clase Graph
    final Graph g1 = new Graph();
    final Graph g2 = new Graph();
    final Graph g3 = new Graph();
    final Graph g4 = new Graph();

    

    //objetos de la clase HashMap para guardar el regreso de los metodos
    //de generación de grafos

    HashMap <String,Double> Nodes = new HashMap <String,Double>();
    HashMap <String,Integer> ErdRny = new HashMap <String,Integer>();
    HashMap <String,Integer> Gilbert = new HashMap <String,Integer>();
    HashMap <String,Integer> SimpleGeo = new HashMap <String,Integer>();
    HashMap <String,Integer> BarabasiAlbert = new HashMap <String,Integer>();

    HashMap <String,Integer> G_BFS = new HashMap <String,Integer>();
    HashMap <String,Integer> G_DFSi = new HashMap <String,Integer>();
    HashMap <String,Integer> G_DFSr = new HashMap <String,Integer>();
    HashMap <String,Double> W_graph = new HashMap <String,Double>();
    Map <String,Double> Dist_nodes = new HashMap<String,Double>();

    Scanner keyboard = new Scanner(System.in);

    int n=500;         //variable numero de nodos
    int m=2000;        //variable numero de aristas
    double p=0.7;     //variable método Gilbert
    int D=25;         //variable método Barabasi Albert
    double r=0.28;    //variable método Geografica Simple
    Boolean d=false;  //variable Dirigido
    Boolean a=false;  //variable Autoconectado
    int v_s=2;        //variable nodo fuente de metodos BFS
    double min=3.0;   //variable mínimo valor aleatorio en aristas
    double max=30.0;  //variable máximo valor aleatroio en aristas
    
    Nodes= g1.node(n);                          //método generador de lista con n nodos     

    ErdRny=g1.genErdosRenyi(n,m,d,a);           
    G_BFS=g1.graphToBFS(ErdRny,Nodes,v_s,d);
    G_DFSi=g1.graphToDFSi(ErdRny,Nodes,v_s,d);
    G_DFSr=g1.graphToDFSr(ErdRny, Nodes, v_s, d);
    W_graph=g1.RandomEdgeValues(ErdRny, min, max);
    Dist_nodes=g1.getDijkstra(W_graph, Nodes, "N4", d);

    try {
	    g1.toViz("ErdRny_n"+n+"_m"+m, d, ErdRny);
  	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
    }/*
    try {
	    g1.toViz("ErdRny_n"+n+"_m"+m+"BFS", d, G_BFS);
  	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
    }
    try {
	    g1.toViz("ErdRny_n"+n+"_m"+m+"DFSi", d, G_DFSi);
  	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
    }
    try {
	    g1.toViz("ErdRny_n"+n+"_m"+m+"DFSr", d, G_DFSr);
  	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
    } */
    try {
	    g1.toVizDoub("ErdRny_n"+n+"_m"+m+"Dijk", d, Dist_nodes);
  	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
  	}

    Nodes= g2.node(n);
    Gilbert=g2.genGilbert(n,m,p,d,a);
    G_BFS=g2.graphToBFS(Gilbert,Nodes,v_s,d);
    G_DFSi=g2.graphToDFSi(Gilbert,Nodes,v_s,d);
    G_DFSr=g2.graphToDFSr(Gilbert, Nodes, v_s, d);
    W_graph=g2.RandomEdgeValues(Gilbert, min, max);
    Dist_nodes=g2.getDijkstra(W_graph, Nodes, "N4", d);

    try {g2.toViz("Gilbert_n"+n+"_p"+p, d, Gilbert);}
    catch(final IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }/*
    try {g2.toViz("Gilbert_n"+n+"_p"+p+"BFS", d, G_BFS);}
    catch(final IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }
    try {g2.toViz("Gilbert_n"+n+"_p"+p+"DFSi", d, G_DFSi);}
    catch(final IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }
    try {
	    g1.toViz("Gilbert_n"+n+"_m"+m+"DFSr", d, G_DFSr);
  	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
    }*/
    try {g2.toVizDoub("Gilbert_n"+n+"_p"+p+"Dijk", d, Dist_nodes);}
    catch(final IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }


    Nodes= g3.node(n);
    SimpleGeo=g3.genSimpleGeo(n,r,d,a);
    G_BFS=g3.graphToBFS(SimpleGeo,Nodes,v_s,d);
    G_DFSi=g3.graphToDFSi(SimpleGeo,Nodes,v_s,d);
    G_DFSr=g3.graphToDFSr(SimpleGeo, Nodes, v_s, d);
    W_graph=g3.RandomEdgeValues(SimpleGeo, min, max);
    Dist_nodes=g3.getDijkstra(W_graph, Nodes, "N4", d);
    
    try {g1.toViz("SimpleGeo_n"+n+"_r"+r, d, SimpleGeo);}
    catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
  	}/*
    try {g1.toViz("SimpleGeo_n"+n+"_r"+r+"BFS", d, G_BFS);}
    catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
    }
    try {g1.toViz("SimpleGeo_n"+n+"_r"+r+"DFSi", d, G_DFSi);}
    catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
    }
    try {
	    g1.toViz("SimpleGeo_n"+n+"_m"+m+"DFSr", d, G_DFSr);
  	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
    }*/
    try {g2.toVizDoub("SimpleGeo_n"+n+"_r"+r+"Dijk", d, Dist_nodes);}
    catch(final IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }

    Nodes= g4.node(n);
    BarabasiAlbert=g4.genBarabasiAlbert(n,D,d,a);
    G_BFS=g4.graphToBFS(BarabasiAlbert,Nodes,v_s,d);
    G_DFSi=g4.graphToDFSi(BarabasiAlbert,Nodes,v_s,d);
    G_DFSr=g4.graphToDFSr(BarabasiAlbert, Nodes, v_s, d);
    W_graph=g4.RandomEdgeValues(BarabasiAlbert, min, max);
    Dist_nodes=g4.getDijkstra(W_graph, Nodes, "N4", d);


    try {g1.toViz("BarabasiAlbert_n"+n+"_d"+D, d, BarabasiAlbert);}
    catch(IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }/*
    try {g1.toViz("BarabasiAlbert_n"+n+"_d"+D+"BFS", d, G_BFS);}
    catch(IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }
    try {g1.toViz("BarabasiAlbert_n"+n+"_d"+D+"DFSi", d, G_DFSi);}
    catch(IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }
    try {
	    g1.toViz("BarabasiAlbert_n"+n+"_m"+m+"DFSr", d, G_DFSr);
  	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
    }*/
    try {g2.toVizDoub("BarabasiAlbert_n"+n+"_D"+D+"Dijk", d, Dist_nodes);}
    catch(final IOException e) {
      System.out.println("An error occurred.");
      e.printStackTrace();
    }

    keyboard.close();
  }
}