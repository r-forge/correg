#' decision tree in a recursive way
#'@export
#'@param data the dataset including the response
#'@param Y the name of the response
#'@param modele vector of the names used for the last tree search (not necessarily choosen)
#'@param kill vector of the names to kill
#'@param index to give a number to the plot
#'@param print boolean to print the tree parameters
#'@param plot boolean to plot the tree
#'@param main the main title if plot=TRUE
#'@param sub the subtitle (if NULL it is automatically added)
#'@param lang the language for the automatic subtitle in the plot
#'
recursive_tree<-function(data=data,Y="Y",modele=NULL,kill=NULL,index=NULL,print=TRUE,plot=TRUE,main=NULL,sub=NULL,lang=c("en","fr")){
   if (is.null(modele)){
      modele=colnames(data)
      modele=modele[modele!=Y]   
   }
   if(!is.null(kill)){
      modele=modele[!modele%in%kill]
   }
   
   formule=as.formula(paste(Y," ~",paste(modele,collapse="+")))
   arbre=rpart(formule,data)
   if(plot){
      opar<-par()
      par( xpd=NA)
      if(is.null(sub)){
         lang=lang[1]
         if(lang=="fr"){
            sub=paste("oui a gauche, non a droite. Moyenne (global", round(arbre$frame$yval[1],digits=3),") et effectif (global ", arbre$frame$n[1],") ")
            vertical= "la hauteur indique la significativite"        
         }else if (lang=="en"){
            sub=paste("True left, False right. Mean (global", round(arbre$frame$yval[1],digits=3),") and effectives (global ", arbre$frame$n[1],") ")
            vertical=" heights indicates significativity"         
         }
      }
      plot(arbre)
      text(arbre, use.n=TRUE)
      title(main = main,ylab=vertical,xlab = sub, col.main = "red", col.lab = gray(.5),cex.main = 1.2, cex.lab = 1.0, font.main = 4, font.lab = 3)
      par( xpd=opar$xpd)
   }
   if(print){
      print(arbre)
   }
   return(list(modele=modele,tree=arbre))
}