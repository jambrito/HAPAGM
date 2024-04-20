hpagm<-function(npart=2,min_by_part=5)
{
#m = matriz com as arestas - colunas 1 e 2 e pesos das arestas a partir da
#coluna 3
#npart = numero de componentes espaciais = grupos com contiguidade
#min_by_part = numero de objetos por grupo
#Md - matriz de distâncias
  
  
  
  
#############################################################################  
escolhe_melhor<-function(sx,Md)
{ Fobj<-rep(0,length(sx))
  for(j in 1:length(sx))
     {clusters<-sx[[j]]
      maxg<-max(clusters)
      grupos<-apply(as.matrix(1:maxg),1,function(cx)  which(clusters==cx))
      Fobj[j]<-sum(sapply(grupos,function(x)  mean(Md[t(comb_n(x,2))])))
   }
return(which.min(Fobj))
}  
  

    
particiona_agm<-function(ga,e,mbp=min_by_part)
{
agm_part<-delete.edges(ga,e)
conexo<-components(agm_part)
if(min(conexo$csize)>=mbp) 
  {cluster<-conexo$membership
  }
else
 {cluster<-NULL
 }
return(cluster)
}  
  
#############################################################################
library(igraph)
library(Rfast)
library(parallel)
library(clusterSim)
nucleos<-detectCores(logical=F)
clust<-makeCluster(nucleos)
clusterEvalQ(clust, library(igraph))
####Leituras#####################################################################
setwd("E:\\Trabs\\Projetos_Finais_ENCE\\Projeto_Nathan_e_Luis_Felipe_2024\\")
library(readxl)
m<-as.matrix(read_excel("distancias_vizinhos_vscc.xlsx"))
y=t(apply(m[,1:2],1,sort))
ix<-which(duplicated.matrix(y)==FALSE)
m<-m[ix,]
Deuclidiana<-as.matrix(read_excel("matriz_euclid_vscc.xlsx"))
DManhattan<-as.matrix(read_excel("matriz_manhat_vscc.xlsx"))
DMahalanobis<-as.matrix(read_excel("matriz_mahal_vscc.xlsx"))
#################################################################################


tempo_cpu<-proc.time()

g<-make_graph(as.vector(t(m[,1:2])),directed = FALSE)
pesos<-m[,3:(3+(ncol(m)-3))]
nnos<-length(union(unique(m[,1]),unique(m[,2])))
if (is.matrix(pesos)==FALSE) {pesos<-as.matrix(pesos)}
AGM<-t(apply(pesos,2,function(w) mst(g,weights=w)))
narestas<-nnos-1
cbe<-comb_n(narestas,npart-1)
cat("Numero de particoes = ",ncol(cbe),"\n")

######################Para Euclidiana#######################
cat("Construindo Partições com a Distância Euclidiana \n")
sol_part=parApply(cl=clust,cbe,2,function(e) particiona_agm(AGM[[1]],e))
npart_viavel<-which(sapply(sol_part,function(x) is.null(x))==FALSE)
npv1<-length(npart_viavel)
solucoes_viaveis<-sol_part[npart_viavel]
if (is.null(solucoes_viaveis)==FALSE)
  {smelhor1<-escolhe_melhor(solucoes_viaveis,Deuclidiana)
   sv<-solucoes_viaveis[[smelhor1]]
  }
 else {sv=rep(0,nnos);smelhor1=0} 
 write.table(sv,
            file=paste("Solucao_Euclidiana_Grupos_VSCC_",npart,"_Cap_",min_by_part,".txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)
######################Para Manhattan#######################
cat("Construindo Partições com a Distância de Manhattan \n")
sol_part=parApply(cl=clust,cbe,2,function(e) particiona_agm(AGM[[2]],e))
npart_viavel<-which(sapply(sol_part,function(x) is.null(x))==FALSE)
npv2<-length(npart_viavel)
solucoes_viaveis<-sol_part[npart_viavel]
if (is.null(solucoes_viaveis)==FALSE)
  {smelhor2<-escolhe_melhor(solucoes_viaveis,DManhattan)
   sv<-solucoes_viaveis[[smelhor2]]
  }
else {sv<-rep(0,nnos);smelhor2=0}
write.table(sv,
            file=paste("Solucao_Manhattan_Grupos_VSCC_",npart,"_Cap_",min_by_part,".txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)


######################Para Mahalanobis#######################
cat("Construindo Partições com a Distância de Mahalanobis \n")
sol_part=parApply(cl=clust,cbe,2,function(e) particiona_agm(AGM[[2]],e))
npart_viavel<-which(sapply(sol_part,function(x) is.null(x))==FALSE)
npv3<-length(npart_viavel)
solucoes_viaveis<-sol_part[npart_viavel]
if (is.null(solucoes_viaveis)==FALSE)
 {smelhor3<-escolhe_melhor(solucoes_viaveis,DMahalanobis)
  sv<-solucoes_viaveis[[smelhor3]]
 }
else {sv=rep(0,nnos);smelhor3=0}
write.table(sv,
            file=paste("Solucao_Mahalanobis_Grupos_VSCC_",npart,"_Cap_",min_by_part,".txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)

tempo_cpu<-(proc.time()-tempo_cpu)[3] 
stopCluster(cl=clust)


return(list(tempo_cpu=tempo_cpu,total_particoes=choose(narestas,npart-1),
            npvs=c(npv1,npv2,npv3),nsm=c(smelhor1,smelhor2,smelhor3)))
}