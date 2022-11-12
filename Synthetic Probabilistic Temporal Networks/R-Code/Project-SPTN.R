library(igraph)
library(centiserve)
library(randomcoloR)
library(igraph,RcolorBrewer)


# ---------------- UTILITY FUNCTIONS -----------------------------------------------------------------------

# Standard Format to display graph
plotGraph <- function(g) {
  #layout_on_sphere
  #layout_with_kk
  set.seed(3014)
  if(vcount(g)!=0) {
    # Plot degree distribution
    plot(degree.distribution(g), log = "xy")
    # Plot graph
    plot(g,layout=layout_with_kk, vertex.label=NA, margin=-0.15,
         vertex.size=10*log(degree(g)+2)*(1/log2(vcount(g)+2)),
         vertex.color="#800000",
         vertex.frame.color="#800000",
         edge.color = "#347E7D",
         #edge.color=ifelse(E(g)$weight==1,"magenta","#347E7D"),
         #edge.width=(1/log2(ecount(g)+2))*ifelse(is.weighted(g),E(g)$weight,1),
         edge.arrow.size=0.5
    )
    print(paste("vertices=",vcount(g)," Edges=",ecount(g)))
  }
}

# Delete dead nodes
deleteNodes <- function(G, t, nLifeFunction) {
  d <- 1
  while (vcount(G)>0 && d<=vcount(G)) {
    v = V(G)[d]
    if(nLifeFunction(G,t,v)<=t-get.vertex.attribute(G, "timestamp",v)){
      G <- delete.vertices(G,c(v))
    }
    d=d+1
  }
  return(G)
}

# Give birth to new nodes
addNodes <- function(G, t, ncrFunction) {
  requiredNodes = vcount(G) + ncrFunction(G,t)
  # If the total required nodes is greater than actual, add rows and columns to adjacency matrix
  if(floor(requiredNodes)>vcount(G)) {
    G <- add_vertices(G,floor(requiredNodes)-vcount(G),"timestamp"=t)
  }
  return(G)
}

# Probabilistic function to add/delete edges between nodes.
addDeleteEdges <- function(G, t, ecf, edf, ewf, selfLoops) {
  # Check if G is not empty
  if(vcount(G)>0){
    #Permutation of vertices
    obj<-ecf(G,t,1,2)
    if(is.list(obj)) {
      G<-obj
    } else{
      for(v1 in V(G)) {
        for(v2 in V(G)) {
          # Check if selfLoop is allowed
          if(selfLoops || v1!=v2){
            p <- ecf(G,t,v1,v2) # Find edge creation probability
            if(!(are.connected(G,v1,v2))){
              if(runif(1)<=p){
                if(is.weighted(G)){
                  G<-add_edges(G,c(v1,v2),"weight"=ewf(G,t,v1,v2))
                  #print(paste("Weights=",E(G)$weight))
                }else{
                  G<-add_edges(G,c(v1,v2))
                }
              }
            }
            if(are.connected(G,v1,v2)){
              p = edf(G,t,v1,v2) # Find edge deletion probability
              aa=ecount(G)
              #print(paste("edges=",aa))
              if(runif(1)<=p) {
                G<- delete_edges(G,(as_ids(E(G)[c(v1) %--% c(v2)])))
              }
            }
          }
        }
      }
    }
  }
  return(G)
}

# Probabilistic function to add/delete edges between nodes.
fixedAddDeleteEdges <- function(G, t, ecf, edf, ewf, selfLoops, fixedEdgeIncrement) {
  # Check if G is not empty
  edgeCount = 0
  if(vcount(G)>0){
    while(edgeCount<fixedEdgeIncrement) {
      # print(paste(edgeCount,"sdfs",fixedEdgeIncrement))
      #Permutation of vertices
      for(v1 in V(G)) {
        for(v2 in V(G)) {
          # Check if selfLoop is allowed
          if(selfLoops || v1!=v2){
            p = ecf(G,t,v1,v2) # Find edge creation probability
            if(!(are.connected(G,v1,v2))){
              if(runif(1)<=p){
                if(is.weighted(G)){
                  if(edgeCount<fixedEdgeIncrement) {
                    G<-add_edges(G,c(v1,v2),"weight"=ewf(G,t,v1,v2))
                    edgeCount=edgeCount+1
                  }
                  #print(paste("Weights=",E(G)$weight))
                }else{
                  if(edgeCount<fixedEdgeIncrement) {
                    G<-add_edges(G,c(v1,v2))
                    edgeCount=edgeCount+1
                  }
                }
              }
            } 
            if(are.connected(G,v1,v2)){
              p = edf(G,t,v1,v2) # Find edge deletion probability
              aa=ecount(G)
              #print(paste("edges=",aa))
              if(runif(1)<=p) {
                G<- delete_edges(G,(as_ids(E(G)[c(v1) %--% c(v2)])))
              }
            }
          }
        }
      }
    }
  }
  return(G)
}


# ----------------------- Temporal Function ---------------------------------------------------------
i <- 0
# Trigger the simulation of temporal graph
temporalNetwork <- function(G, ncrFunction, nLifeFunction, ecf, edf, ewf, t_start, t_end, unit_time, selfLoops, fixedEdgeIncrement) {
  # Generating Empty MetaMap
  metamap<-list(c(),c(),c(),c(),c(),c(),c(),c(),c(),c())
  names(metamap)<-c("Time", "Nodes","Edges","AvgDegree","MaxDegree","Components","AvgPathLength","Diameter","ALCC","GCC")
  # initializing timestamp for each vertex
  for(v in V(G)) {
    G <- set_vertex_attr(G, "timestamp",v, t_start)
  }
  plotGraph(G) # Plot of initial graph
  # For time t = start to upperbound
  for(t in t_start+1:t_end){
    print(paste("t=",t))
    Sys.sleep(unit_time)
    #G <- deleteNodes(G, t, nLifeFunction)
    G <- addNodes(G, t, ncrFunction)
    if(fixedEdgeIncrement==-1) {
      G <- addDeleteEdges(G, t, ecf, edf, ewf, selfLoops)
    } else {
      G <- fixedAddDeleteEdges(G, t, ecf, edf, ewf, selfLoops, fixedEdgeIncrement)
    }
    if(t%%1==0) {
      metamap <- metaDataGenerator(G, metamap, t)
      plotGraph(G)
    }
  }
  return(metamap)
}

#---------------------------------------INPUTS-----------------------------------

#------------------Example Network-----------------------------------------------
# # Node Creation Rate function passed as an input
# # Return a node creation rate per unit time t
# ncrFunction <- function(G, t) {
#   return (exp(t/15)) # Creation exponentially depends on time.
# }
# # Node Life function passed as an input
# # Return the lifetime of a node in standard unit time scale (t)
# nLifeFunction <- function(G,t,v) {
#   birthTime <- get.vertex.attribute(G, "timestamp",v)
#   return (100*runif(1,0.7,1)) # Age anywhere between 70-100
# }
# # Edge Creation function passed as an input
# # Return a probability dependent on any of the graph/node attributes
# edgeCreationFunction <- function(G,t,v1,v2) {
#   # Inversely proportional to age, degree of the vertices.
#   birthTime1<-get.vertex.attribute(G, "timestamp",v1)
#   birthTime2<-get.vertex.attribute(G, "timestamp",v2)
#   p1<-1/(1+exp(-(t-birthTime1+5)))
#   p2<-1/(1+exp(-(t-birthTime2+5)))
#   p<-(p1*p2)/2
#   d<-1/(degree(G,v1,mode="all")*degree(G,v2,mode="all"))
#   return(p*d/vcount(G))
# }
# # Edge Deletion function passed as an input
# # Return a probability dependent on any of the graph/node attributes
# edgeDeletionFunction <- function(G,t,v1,v2) {
#   return(1/(1+log(t+1))) # Inversely proportional ot log(t)
# }

#------------------Erdos Renyi----------------------
# ncrFunction <- function(G, t) {
#   return (150)
# }
# nLifeFunction <- function(G,t,v) {
#  
#   return (2000)
# }
# edgeCreationFunction <- function(G,t,v1,v2) {
#   if(get.vertex.attribute(G, "timestamp", v1)==t) {
#     return(0.1)
#   } else {
#     return(0)
#   }
# }
# edgeDeletionFunction <- function(G,t,v1,v2) {
#   return(0)
# }


#------------------Barabasi Albert------------------
# ncrFunction <- function(G, t) {
#   return(1) # Barabasi has 1 node creation per unit time
# }
# nLifeFunction <- function(G,t,v) {
#   return(2000) # no node deletion. Infinite life.
# }
# edgeCreationFunction <- function(G,t,v1,v2) {
#   birthTime<-get.vertex.attribute(G, "timestamp",v1)
#   p<-1 # Power
#   if(birthTime==t){
#     if(degree(G, v2, mode = "in")!=0){
#     z<-(degree(G, v2, mode = "in")*p+1)/((sum(degree(G, mode = "in"))*p)+vcount(G))
#     return(z) # probability of edge creation
#     } else {
#       return(1) # zero.appeal is used if indegree=0
#     }
#   } else {
#     return(0) # for old nodes no edge creation
#   }
# }
# edgeDeletionFunction <- function(G,t,v1,v2) {
#   return(0)
# }


#------------------R-MAT----------------------------
ncrFunction <- function(G, t) {
  return (100) # 100 nodes per unit time
}
nLifeFunction <- function(G,t,v) {
  return(20000) # Infinite life
}
edgeCreationFunction <- function(G,t,v1,v2) {
  # This function has the R-MAT logic
  B<-RMATedgeCreationFunction(G, v1, v2, t)
  return(B) # Returns graph object
}
edgeDeletionFunction <- function(G,t,v1,v2) {
  return(0) # No edge deletion
}
#
RMATedgeCreationFunction <- function(G,t,v1,v2) {
  prob<-c(0.4,0.6,0.8,1)
  start<-1
  end<-vcount(G)
  sec<-c(1,1,vcount(G),vcount(G))
  section = 0
  ed<-ecount(G)
  while(ecount(G)-ed!=100) {
    # print(paste("i",i))
    G<-recursive(prob, sec, get.adjacency(G))
  }
  return(G)

}
recursive<-function(prob, sec, A) {
  if(sec[1]==sec[3] && sec[2]==sec[4]) {
    if(sec[1]!=sec[2]) {
      G<-add_edges(graph_from_adjacency_matrix(A),c(sec[1],sec[2]))
    } else {
      G<-graph_from_adjacency_matrix(A)
    }
    return(G)
  } else {
    p<-runif(1)
    for(i in 1:length(prob)) {
      if(p<=prob[i]) {
        section = i
        break
      }
    }
      if(section==1) {
        sec<-c(sec[1],sec[2],floor((sec[3]+sec[1])/2),floor((sec[4]+sec[2])/2))
      } else if(section==2) {
        sec<-c(sec[1],floor((sec[2]+sec[4])/2),floor((sec[3]+sec[1])/2),sec[4])
      } else if(section==3) {
        sec<-c(floor((sec[1]+sec[3])/2),sec[2],sec[3],floor((sec[2]+sec[4])/2))
      } else if(section==4) {
        sec<-c(floor((sec[1]+sec[3])/2),floor((sec[3]+sec[1])/2),sec[3],sec[4])
      }
    recursive(prob, sec, A)
  }
}


#------------------Planted Cliques------------------
# ncrFunction <- function(G, t) {
#   return (2) # 2 nodes per unit time
# }
# nLifeFunction <- function(G,t,v) {
#   return(2000) # Large or infinite life
# }
# edgeCreationFunction <- function(G,t,v1,v2) {
#   v = c()
#   if(t%%20==0) {
#     v = sample.int(vcount(G), 4)
#     #v = runif(10, 1, vcount(G))
#     A = get.adjacency(G)
#     for(i in v) {
#       for(j in v) {
#         if(i!=j) {
#           A[i,j]=1
#           A[j,i]=1
#           E(G)$color="black"
#         }
#       }
#     }
#     g<-graph_from_adjacency_matrix(A)
#     V(g)$color = "#800000"
#     for(vv in v){
#       V(g)$color[vv]="blue"
#     }
#     
#   }
#  
#   if(get.vertex.attribute(G, "timestamp",v1)==t) {
#       return(0.1)
#   } else {
#       return(0)
#     }
#     
# }
# edgeDeletionFunction <- function(G,t,v1,v2) {
#   return(0)
# }


#------------------Signed Networks------------------
# ncrFunction <- function(G, t) {
#   return (2) # Two nodes per unit time
# }
# nLifeFunction <- function(G,t,v) {
#   return (5) # Node life = 5 units time for each node
# }
# edgeCreationFunction <- function(G,t,v1,v2) {
#   return(0.1) # Edge creation probability = -1
# }
# edgeDeletionFunction <- function(G,t,v1,v2) {
#   return(0) # Edge deletion probability = -1
# }
# edgeWeightFunction <- function(G,t,v1,v2) {
#   if(runif(1)<0.5) {
#     return(1) # Positive Edge
#   } else {
#     return(2) # Negative Edge
#   }
# }

#------------------------------------MAIN FUNCTION----------------------------------
# Main function - Take Inputs
mainFunction <- function() {
  t_start<-0 # Start time
  t_end<-10 # End time
  A <- matrix(c(0,0,0,0),nrow=2) # Adjacency Matrix of Initial graph. First row is the time stamp of the nodes.
  directionType = "directed"
  G <- make_empty_graph(n = 2, directed = TRUE)
  #G <- graph_from_adjacency_matrix(A, mode=directionType, weighted = NULL)
  plot(G)
  selfLoops = FALSE
  unit_time<-0 # Unit time = 1 sec
  
  # Trigger the temporal graph generation function
  metamap <- temporalNetwork(G ,ncrFunction, nLifeFunction, edgeCreationFunction, edgeDeletionFunction, edgeWeightFunction, t_start, t_end, unit_time, selfLoops,-1)
  print(paste("MetaMap = ", metamap))
  
  # ba<-barabasi.game(500, 1, m=2, zero.appeal = 1, directed = TRUE, start.graph = G)
  # plotGraph(ba)
  # print(paste(ecount(ba),"bara-ver",vcount(ba)))
  # print(paste("edges=",ecount(G),"degree",degree(G)))
  # er <- erdos.renyi.game(1000, 0.1, loops = FALSE)
  # plotGraph(er)
  # metamap <- metaDataGenerator(er,metamap,t)
  
}

metaDataGenerator<-function(G, metamap, t) {
  # c("Time", "Nodes","Edges","AvgDegree","MaxDegree","Components","AvgPathLength","Diameter","ALCC","GCC")
  metamap$Time <- c(metamap$Time,t)
  metamap$Nodes <- c(metamap$Nodes,vcount(G))
  metamap$Edges <- c(metamap$Edges, ecount(G))
  metamap$AvgDegree <- c(metamap$AvgDegree,mean(degree(G, mode = "all", loops = FALSE, normalized = FALSE)))
  metamap$MaxDegree <- c(metamap$MaxDegree,max(degree(G, mode = "all", loops = FALSE, normalized = FALSE)))
  metamap$Components <- c(metamap$Components, count_components(G, mode ="weak"))
  metamap$AvgPathLength <- c(metamap$AvgPathLength, mean_distance(G, directed = is.directed(G), unconnected = !is.connected(G, mode = "weak")))
  metamap$Diameter <- c(metamap$Diameter, diameter(G, directed = is.directed(G), unconnected = !is.connected(G, mode = "weak"), weights = NULL))
  metamap$ALCC <- c(metamap$ALCC, transitivity(G, type = "average"))
  metamap$GCC <- c(metamap$GCC, transitivity(G, type = "global") )
  return(metamap)
}

plotLineGraph<-function(x,xlabel,y,ylabel){
  plot(x, y, type = "o", frame = FALSE, pch = runif(1,1,15), 
       col = runif(1,1,15), xlab =xlabel , ylab = ylabel)
}

# Starting Point
mainFunction()






