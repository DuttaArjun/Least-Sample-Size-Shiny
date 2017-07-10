#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  #Header or Title Panel
  titlePanel(title = h4("Investigation of Convergence to Normality of Quantile Measures From Continuous Distributions with R Shiny", align = "center")),
  #Sidebar Layout
  sidebarLayout(
    #Sidebar Panel
    sidebarPanel(
      selectInput("Dist", "Select the Distribution", choice = c("Exponential (rate = m)" = 1, "LogNormal" = 2, "Cauchy" = 3, "Normal" = 4, "Beta" = 5, "Uniform" = 6), selected = "Exponential"),
      selectInput("QM", "Select the Quantile Measure", choice = c("Median" = 1, "Quartile Deviation" = 2, "Coefficient of Quartile Deviation" = 3, "Bowley Measure of Skewness" = 4, "Percentile Measure of Kurtosis" = 5), selected = "Median"),
      conditionalPanel(  condition = "input.Dist == '1'",
                         numericInput("k", "State the parameter m of Exponential Distribution",value = 1)),
      conditionalPanel(  condition = "input.Dist == '5'",
                         numericInput("a", "State the parameter m of Beta Distribution",value = 1),
                         numericInput("b", "State the parameter n of Beta Distribution",value = 1)),
      sliderInput("alpha", "Select the level of Significance ", min = 0, max = 1,step = 0.05, value = 0.05),               
      numericInput("r", "State the number of sample",value = 1000),
      br(),
      br(),
      h4("Powered By "),img(src = 'R.png', height='50px',width='50px')
    ),
    #Main Panel
    mainPanel(h4("Graphs for Different Quantile Measures", align = "center"), 
              plotOutput("plot"),
              textOutput("dnc")
              
    )
  )
)
# Define server logic required to draw a histogram
server <- function(input, output) {
  output$dnc = renderText({
    if ((input$Dist == "4" && input$QM == "3"))
      print("Coefficient of Quantile Deviation does not exist for Normal Distribution")
    else if ((input$Dist == "3" && input$QM == "3"))
      print ("Coefficient of Quantile Deviation does not exist for Cauchy Distribution")
  })
  # Median
  output$plot = renderPlot({
    if (input$QM == "1"){
      n = 2  
      pv = 0
      while(pv<input$alpha)
      {
        med=vector()
        switch(input$Dist, "1" = (for(i in 1:input$r) med[i]=median(rexp(n,input$k))),
               "2" = (for(i in 1:input$r) med[i]=median(rlnorm(n))),
               "3" = (for(i in 1:input$r) med[i]=median(rcauchy(n))), 
               "4" = (for(i in 1:input$r) med[i]=median(rnorm(n))), 
               "5" = (for(i in 1:input$r) med[i]=median(rbeta(n,input$a,input$b))), 
               "6" = (for(i in 1:input$r) med[i]=median(runif(n)))
        )
        pv=as.numeric(shapiro.test(med)[2])
        n=n+1
      }
      hist((med-mean(med))/sd(med),freq=F,main=paste("Median,n=",n-1))
      curve(dnorm,add=T)
    }
    #Quantile Deviation
    else if (input$QM == "2"){
      n = 2  
      pv = 0
      while(pv<input$alpha)
      {
        qd=vector()
        switch(input$Dist, "1" = (for(i in 1:input$r) qd[i]=(quantile(rexp(n,input$k),0.75)-quantile(rexp(n,input$k),0.25))/2), 
               "2" = (for(i in 1:input$r) qd[i]=(quantile(rlnorm(n),0.75)-quantile(rlnorm(n),0.25))/2), 
               "3" = (for(i in 1:input$r) qd[i]=(quantile(rcauchy(n),0.75)-quantile(rcauchy(n),0.25))/2),
               "4" = (for(i in 1:input$r) qd[i]=(quantile(rnorm(n),0.75)-quantile(rnorm(n),0.25))/2),
               "5" = (for(i in 1:input$r) qd[i]=(quantile(rbeta(n,input$a,input$b),0.75)-quantile(rbeta(n,input$a,input$b),0.25))/2),
               "6" = (for(i in 1:input$r) qd[i]=(quantile(runif(n),0.75)-quantile(runif(n),0.25))/2)
        )
        pv=as.numeric(shapiro.test(qd)[2])
        n=n+1
      }
      hist((qd-mean(qd))/sd(qd),freq=F,main=paste("QD,n=",n))
      curve(dnorm,add=T)
    }
    #Coeffcient of Quantile Deviation
    else if (input$QM == "3"){
      n = 2  
      pv = 0
      while(pv<input$alpha)
      {
        cqd=vector()
        switch(input$Dist, "1" = (for(i in 1:input$r) cqd[i]=as.numeric((quantile(rexp(n,input$k),0.75)-quantile(rexp(n,input$k),0.25))/2)/median(rexp(n,input$k))), 
               "2" = (for(i in 1:input$r) cqd[i]=as.numeric((quantile(rlnorm(n),0.75)-quantile(rlnorm(n),0.25))/2)/median(rlnorm(n))),
               "5" = (for(i in 1:input$r) cqd[i]=as.numeric((quantile(rbeta(n,input$a,input$b),0.75)-quantile(rbeta(n,input$a,input$b),0.25))/2)/median(rbeta(n,input$a,input$b))),
               "6" = (for(i in 1:input$r) cqd[i]=as.numeric((quantile(runif(n),0.75)-quantile(runif(n),0.25))/2)/median(runif(n)))
        )
        pv=as.numeric(shapiro.test(cqd)[2])
        n=n+1
      }
      hist((cqd-mean(cqd))/sd(cqd),freq=F,main=paste("CQD,n=",n))
      curve(dnorm,add=T)
    }
    #Bowleys measure of Skewness
    else if (input$QM == "4")
    {
      n = 2
      pv = 0
      while(pv<input$alpha)
      {
        sk=vector()
        switch(input$Dist, 
               "1" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(rexp(n, input$k),0.75)-quantile(rexp(n,input$k),0.25))/2)
                 sk[i]=as.numeric((quantile(rexp(n, input$k),0.75)+quantile(rexp(n,input$k),0.25)-2*median(rexp(n,input$k)))/(2*qd))
               }                      ),
               "2" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(rlnorm(n),0.75)-quantile(rlnorm(n),0.25))/2)
                 sk[i]=as.numeric((quantile(rlnorm(n),0.75)+quantile(rlnorm(n),0.25)-2*median(rlnorm(n)))/(2*qd))
               }                      ),
               "3" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(rcauchy(n),0.75)-quantile(rcauchy(n),0.25))/2)
                 sk[i]=as.numeric((quantile(rcauchy(n),0.75)+quantile(rcauchy(n),0.25)-2*median(rcauchy(n)))/(2*qd))
               }                      ),
               "4" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(rnorm(n),0.75)-quantile(rnorm(n),0.25))/2)
                 sk[i]=as.numeric((quantile(rnorm(n),0.75)+quantile(rnorm(n),0.25)-2*median(rnorm(n)))/(2*qd))
               }                      ),
               "5" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(rbeta(n,input$a,input$b),0.75)-quantile(rbeta(n,input$a,input$b),0.25))/2)
                 sk[i]=as.numeric((quantile(rbeta(n,input$a,input$b),0.75)+quantile(rbeta(n,input$a,input$b),0.25)-2*median(rbeta(n,input$a,input$b)))/(2*qd))
               }                      ),
               "6" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(runif(n),0.75)-quantile(runif(n),0.25))/2)
                 sk[i]=as.numeric((quantile(runif(n),0.75)+quantile(runif(n),0.25)-2*median(runif(n)))/(2*qd))
               }     )
        )
        pv=as.numeric(shapiro.test(sk)[2])
        n = n+1
      }
      hist((sk-mean(sk))/sd(sk),freq=F,main=paste("SK,n=",n))
      curve(dnorm,add=T)
    }
    #Percentile Measure of Skewness
    else if(input$QM == "5")
    {
      n = 2
      pv = 0
      while(pv<input$alpha)
      {
        kp=vector()
        switch(input$Dist, 
               "1" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(rexp(n, input$k),0.75)-quantile(rexp(n,input$k),0.25))/2)
                 kp[i]=as.numeric(qd/(quantile(rexp(n, input$k),0.90)-quantile(rexp(n, input$k),0.10)))
               }                      ),
               "2" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(rlnorm(n),0.75)-quantile(rlnorm(n),0.25))/2)
                 kp[i]=as.numeric(qd/(quantile(rlnorm(n),0.90)-quantile(rlnorm(n),0.10)))
               }                      ),
               "3" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(rcauchy(n),0.75)-quantile(rcauchy(n),0.25))/2)
                 kp[i]=as.numeric(qd/(quantile(rcauchy(n),0.90)-quantile(rcauchy(n),0.10)))
               }                      ),
               "4" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(rnorm(n),0.75)-quantile(rnorm(n),0.25))/2)
                 kp[i]=as.numeric(qd/(quantile(rnorm(n),0.90)-quantile(rnorm(n),0.10)))
               }                      ),
               "5" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(rbeta(n,input$a,input$b),0.75)-quantile(rbeta(n,input$a,input$b),0.25))/2)
                 kp[i]=as.numeric(qd/(quantile(rbeta(n,input$a,input$b),0.90)-quantile(rbeta(n,input$a,input$b),0.10)))
               }                      ),
               "6" = (for(i in 1:input$r) 
               {
                 qd=as.numeric((quantile(runif(n),0.75)-quantile(runif(n),0.25))/2)
                 kp[i]=as.numeric(qd/(quantile(runif(n),0.90)-quantile(runif(n),0.10)))
               }                      )
        )
        pv=as.numeric(shapiro.test(kp)[2])
        n = n+1
      }
      hist((kp-mean(kp))/sd(kp),freq=F,main=paste("KP,n=",n))
      curve(dnorm,add=T)
    }
  })
}
# Run the application 
shinyApp(ui = ui, server = server)