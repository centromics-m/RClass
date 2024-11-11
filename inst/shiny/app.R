shinyApp(
ui<-fluidPage(sliderInput("iter", label="simCLT", min = 0, max = 1000, value = 50), plotOutput("plot")),

server<-function(input, output){ 
  my.clt.sim <- function(rd, iter=1000, n=10) {
    par(mfrow=c(1,2))
    for(i in names(rd)) {
      hist( rd2<- replicate(iter, mean(sample( rd[[i]] , n))  ), main= paste0(i, '\n iter = ', iter, ', n = ', n ), xlab='',  probability = T  ) 
      lines(density(rd2), col=2)
    }
    par(mfrow=c(1,1))
  }
  
  rd <- list(
    runif(100, 1, 100),
    rchisq(100, 3),
    rbinom(100, 10, 1/6), 
    rhyper(100, 10, 90, 10) 
  )   
  names(rd) <- c('Uniform', 'Chi-square', 'Binominal', 'Hypergeometric')  
  
  v<-reactiveValues();
  observeEvent(input$iter,    { v$data<-input$iter } )
  output$plot <- renderPlot({ 
    if(is.null(v$data)) return()
    
    my.clt.sim(rd, v$data)
    
    
  }) }
)
