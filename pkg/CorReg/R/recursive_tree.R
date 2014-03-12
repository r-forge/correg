#' decision tree in a recursive way
#'@export
#'@param data the dataset including the response
#'@param Y the name of the response
#'@param modele vector of the names used for the last tree search (not necessarily choosen)
#'@param kill vector of the names to kill
#'@param index to give a number to the plot
#'@param print boolean to print the tree parameters
#'@param plot boolean to plot the tree
#'
recursive_tree<-function(data=data,Y="Y",modele=NULL,kill=NULL,index=NULL,print=TRUE,plot=TRUE){
   if (is.null(modele)){
      modele=colnames(data)
      modele=modele[modele!=Y]   
   }
   if(!is.null(kill)){
      modele=modele[!modele%in%kill]
   }
   
   formule=as.formula(paste(Y," ~",paste(modele,collapse="+")))
   arbre=rpart(formule,data);par( xpd=NA)
   if(plot){
      sub=paste("oui a gauche, non a droite. Moyenne (global", round(arbre$frame$yval[1],digits=3),") et effectif (global ", arbre$frame$n[1],") \n la hauteur indique la significativite")
      plot(arbre,sub=sub)
      text(arbre, use.n=TRUE)
      #title(main = titre,ylab=vertical,xlab = paste(infos,"\n Oui Ã  gauche, non Ã  droite.",reponse), col.main = "red", col.lab = gray(.5),cex.main = 1.2, cex.lab = 1.0, font.main = 4, font.lab = 3
   }
   if(print){
      print(arbre)
   }
   return(list(tree=arbre,modele=modele))
}