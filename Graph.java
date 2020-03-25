import java.util.*;
import java.io.*;
import java.nio.file.DirectoryIteratorException;

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
        if (dir) {
          if (Auto && (N_i == N_f)) {
            G_ErdosRenyi.put("N" + N_i + "->N" + N_f, i);
            System.out.println(N_i + "->" + N_f);
          }
          if (!Auto && N_i != N_f) {
            G_ErdosRenyi.put("N" + N_i + "->N" + N_f, i);
            System.out.println(N_i + "->" + N_f);
          }
        }
        if (!dir && !G_ErdosRenyi.containsKey("N" + N_f + "--N" + N_i)) {
          if (Auto && (N_i == N_f)) {
            G_ErdosRenyi.put("N" + N_i + "--N" + N_f , i);
            System.out.println(N_i + "--" + N_f);
          }
          if (!Auto && N_i != N_f) {
            G_ErdosRenyi.put("N" + N_i + "--N" + N_f , i);
            System.out.println(N_i + "--" + N_f);
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

  public static HashMap<String, Integer> genGilbert(final int n, final int m, final double p, final Boolean d,
      final Boolean a) {
    final HashMap<String, Integer> G_Gilbert = new HashMap<String, Integer>();
    HashMap<String, Double> Nodos_glbr = new HashMap<String, Double>();
    Nodos_glbr = node(n);
    // variables locales
    final int num_nod = Nodos_glbr.size();
    final int edgs = m;
    final double prbl = p;
    final Boolean dir = d;
    final Boolean Auto = a;
    int edgs_cnt = 0;
    // ciclos para recorrer todas las parejas posibles de nodos y asignar arista
    // si numero aleatorio es superior a p
    for (int i = 0; i < num_nod; i++) {
      for (int j = 0; j < num_nod; j++) {
        if (G_Gilbert.containsKey("N" + i + "->N" + j)) { } 
        else {
          final double rnd = Math.random();
          if (rnd > prbl) {
            if (dir) {
              if (Auto && i == j) {
                G_Gilbert.put("N" + i + "->N" + j, edgs_cnt);
              }
              if (!Auto && i != j) {
                G_Gilbert.put("N" + i + "->N" + j, edgs_cnt);
              }
            }
            if (!dir && !G_Gilbert.containsKey("N" + j + "--N" + i )) {
              if (Auto && i == j) {
                G_Gilbert.put("N" + i + "--N" + j , edgs_cnt);
              }
              if (!Auto && i != j) {
                G_Gilbert.put("N" + i + "--N" + j , edgs_cnt);
              }
            }
          }
          edgs_cnt = +edgs_cnt;
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
  public static HashMap<String, Integer> genSimpleGeo(final int n, final double r, final Boolean d, final Boolean a) {
    final HashMap<String, Integer> G_SimpleGeo = new HashMap<String, Integer>();
    HashMap<String, Double> N_x = new HashMap<String, Double>();
    HashMap<String, Double> N_y = new HashMap<String, Double>();

    final int N = n;
    int edgs_cnt;
    final double B = r;
    double R = 0;
    final Boolean Auto = a;
    final Boolean dir = d;
    N_x = node(N);
    N_y = node(N);

    double x1 = 0;
    double x2 = 0;
    final double y1 = 0;
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
          if (G_SimpleGeo.containsKey("N" + i + "->N" + j)) {
          } else {
            if (dir) {
              if (Auto && i == j) {
                G_SimpleGeo.put("N" + i + "->N" + j, edgs_cnt);
              }
              if (!Auto && i != j) {
                G_SimpleGeo.put("N" + i + "->N" + j, edgs_cnt);
              }
            }
            if (!dir && !G_SimpleGeo.containsKey("N" + j + "--N" + i)) {
              if (Auto && i == j) {
                G_SimpleGeo.put("N" + i + "--N" + j , edgs_cnt);
              }
              if (!Auto && i != j) {
                G_SimpleGeo.put("N" + i + "--N" + j , edgs_cnt);
              }
            }
          }
        }
        edgs_cnt = edgs_cnt + 1;
      }
    }

    return G_SimpleGeo;
  }
//arreglar la duplicidad en no dirigidos
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
            if (dir) {
              if (Auto && ("N" + i == "N" + j)) {
                G_BarabasiAlbert.put("N" + i + "->N" + j, edgs_cnt);
                double temp_deg = nodos.get("N" + i);
                nodos.put("N" + i, temp_deg + 1);
                temp_deg = nodos.get("N" + j);
                nodos.put("N" + j, temp_deg + 1);
                edgs_cnt = edgs_cnt + 1;
              }
              if (!Auto && ("N" + i != "N" + j)) {
                G_BarabasiAlbert.put("N" + i + "->N" + j, edgs_cnt);
                double temp_deg = nodos.get("N" + i);
                nodos.put("N" + i, temp_deg + 1);
                temp_deg = nodos.get("N" + j);
                nodos.put("N" + j, temp_deg + 1);
                edgs_cnt = edgs_cnt + 1;
              }
            }
            if (!dir && !G_BarabasiAlbert.containsKey("N" + j + "--N" + i)) {
              if (Auto && ("N" + i == "N" + j)) {
                G_BarabasiAlbert.put("N" + i + "--N" + j , edgs_cnt);
                double temp_deg = nodos.get("N" + i);
                nodos.put("N" + i, temp_deg + 1);
                temp_deg = nodos.get("N" + j);
                nodos.put("N" + j, temp_deg + 1);
                edgs_cnt = edgs_cnt + 1;
              }
              if (!Auto && ("N" + i != "N" + j)) {
                G_BarabasiAlbert.put("N" + i + "--N" + j , edgs_cnt);
                double temp_deg = nodos.get("N" + i);    
                nodos.put("N" + i, temp_deg + 1);
                temp_deg = nodos.get("N" + j);
                nodos.put("N" + j, temp_deg + 1);
                edgs_cnt = edgs_cnt + 1;
              }
            }
          }
        }
      }

    }

    return G_BarabasiAlbert;
  }

  // metodo para guardar la lista generada por los diferentes metodos en archivo
  // .viz
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
        writer.write("\t" + k + ";\n");
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
          writer.write("\t" + k + ";\n");
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

//metodo para separar las llaves de las aristas de cada metodo
//preservamos el primer nodo
  public String unShambleFirst(String bloody_key){
    String str = bloody_key;
    String[] blStr = str.split("->", 2); 
   
    return blStr[0];
  }

  public String unShambleSecond(String bloody_key, Boolean Dir){
    Boolean Dirigido=Dir;
    String str = bloody_key;
    String[] blStr = str.split("->", 2);

    return blStr[1];
  }

//metodo BFS de forma iterativa.

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
            G_BFS.put(key1+"--"+key2, edges_counter);
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
            G_BFS.put(key1+"->"+key2, edges_counter);
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
          G_DFSi.put(key2+"->"+key1, 3);
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
          G_DFSi.put(key2+"--"+key1, 3);
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
  Map<String,Integer> map_to_hold_the_bloody_and_uniterable_HashMap_edges_plus_plus = new HashMap<String,Integer>(Graph_in);
  Map<String,Double> another_map_to_hold_the_other_HashMap_nodes = new HashMap<String,Double>(Nodos_G_in);
  Map<String,Boolean> NotDiscover_map = new HashMap<String,Boolean>();
  Stack<String> stack = new Stack<String>();

  int Nodo_previo=N_source;
  int edges_counter=0;
  int Numero_Nodos_descubiertos=1;
  int nodos_restantes;
  String key1;
  Boolean Dirigido=Dir;

  //quiero recorrer la lista de nodos del grafo para ir buscando desde el nodo fuente u
  //si existen las aristas (u,v) en el grafo e ir añadiendo los nodos al stack
  //una vez descubiertos todos los nodos conectados con u llamar nuevamente a la función  
  for(String key2: another_map_to_hold_the_other_HashMap_nodes.keySet()){
    map_to_hold_the_bloody_and_uniterable_HashMap_edges_plus_plus.put(key2, 1);
  }
  //pero la función debe devolver la lista de nodos descubiertos el stack para el orden
  // de recorrido en el grafo y la lista de aristas del grafo DFS
  // no se como return puede devolver mas de un elemento (HashMap, Stack, map y nuevo nodo previo )
  // para poder realizar la recursividad.
  key1="N"+Nodo_previo;
  map_to_hold_the_bloody_and_uniterable_HashMap_edges_plus_plus.put(key1, 0);
  
  

  return G_DFSr;
}
//Tambien intenté crear una clase para que los objetos tengan los 4 atributos que necesito
//y con un metodo poder modificar éstos pero  
static Map<String, Boolean> Small_DFS_recursive(Map<String,Integer> Map_edges, String v){
  //HashMap para manipular el HashMap entrante dentro de la función recursiva
  Map<String,Integer> map_to_hold_the_bloody_and_uniterable_HashMap_edges_PP = new HashMap<String,Integer>();
  
  map_to_hold_the_bloody_and_uniterable_HashMap_edges_PP=Map_edges;
  
  String key1=v;

 
  for(String key2: map_to_hold_the_bloody_and_uniterable_HashMap_edges_PP.keySet()){

    if(map_to_hold_the_bloody_and_uniterable_HashMap_edges_PP.containsKey(key1+"--"+key2)){
      //si llamo a esta misma función debo entregarle el HashMap con las aristas, la lista con los 
      //nodos visitados, otro HashMap con el DFS hasta el momento y el nodo fuente. Pero no se como 
      //sacarlos de esta función para que loes entregue a si misma
    }
  }

  return Small_DFS_recursive(map_to_hold_the_bloody_and_uniterable_HashMap_edges_PP, key1);
}


  public static void main(final String[] args){
    // creacion de objetos de la clase Graph
    final Graph g1 = new Graph();
    final Graph g2 = new Graph();
    final Graph g3 = new Graph();
    final Graph g4 = new Graph();

    //variable numero de nodos
    int n=0;
    //variable numero de aristas
    int m;
    //variable metodo Gilbert
    double p;
    //variable metodo Geografica Simple
    double r=0.08;
    //variable metodo Barabasi Albert
    int D=5;
    //variable Dirigido
    Boolean d=false;
    //variable Autoconectado
    Boolean a=false;
    //variable para el nodo fuente de metodos BFS
    int v_s=0;
    //objetos de la clase HashMap para guardar el regreso de los metodos
    HashMap <String,Double> Nodes = new HashMap <String,Double>();
    HashMap <String,Integer> ErdRny = new HashMap <String,Integer>();
    HashMap <String,Integer> Gilbert = new HashMap <String,Integer>();
    final HashMap <String,Integer> SimpleGeo = new HashMap <String,Integer>();
    final HashMap <String,Integer> BarabasiAlbert = new HashMap <String,Integer>();
    HashMap <String,Integer> G_BFS = new HashMap <String,Integer>();
    HashMap <String,Integer> G_DFSi = new HashMap <String,Integer>();
    HashMap <String,Integer> G_DFSr = new HashMap <String,Integer>();

    final Scanner keyboard = new Scanner(System.in);
/*
    System.out.println("Numero de nodos");
    n=keyboard.nextInt();
/*
    System.out.println("Numero de aristas ");
    m=keyboard.nextInt();

    System.out.println("Probabilidad para metodo Gilbert");
    p=keyboard.nextDouble();
/*
    System.out.println("Radio de la vecindad para metodo Geografico simple");
    r=keyboard.nextDouble();

    System.out.println("Parametro d para metodo Barabasi Albert");
    D=keyboard.nextInt();

    System.out.println("Dirigido");
    d=keyboard.nextBoolean();

    System.out.println("Autoconectado");
    a=keyboard.nextBoolean();
/
    System.out.println("nodo fuente ");
    v_s=keyboard.nextInt();
*/  n=100;
    m=200;
    p=0.8;
    v_s=2;
    
    Nodes= g1.node(n);
    

    //ErdRny=g1.genErdosRenyi(n,m,d,a);
    Gilbert=g2.genGilbert(n,m,p,d,a);
    G_BFS=g2.graphToBFS(Gilbert,Nodes,v_s,d);
    G_DFSi=g2.graphToDFSi(Gilbert,Nodes,v_s,d);
   /* 
    SimpleGeo=g3.genSimpleGeo(n,r,d,a);
    BarabasiAlbert=g4.genBarabasiAlbert(n,D,d,a);

	try {
	    g1.toViz("ErdRny_n"+n+"_m"+m, d, ErdRny);
	}catch(IOException e) {
	    System.out.println("An error occurred.");
	    e.printStackTrace();
	}
*/
	try {
	    g1.toViz("Gilbert_n"+n+"_p"+p, d, Gilbert);
	}catch(final IOException e) {
	    System.out.println("An error occurred.");
      e.printStackTrace();
  }
  
  try {
    g1.toViz("Gilbert_n"+n+"_p"+p+"BFS", d, G_BFS);
  }catch(final IOException e) {
    System.out.println("An error occurred.");
    e.printStackTrace();
  }
  
  try {
  g1.toViz("Gilbert_n"+n+"_p"+p+"DFSi", d, G_DFSi);
  }catch(final IOException e) {
  System.out.println("An error occurred.");
  e.printStackTrace();
  }
/*
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
*/
    keyboard.close();
  }
}