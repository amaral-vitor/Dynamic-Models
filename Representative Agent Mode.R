#install.packages("shiny")
#install.packages("ggplot2")
#install.packages("deSolve")

library(shiny)
library(ggplot2)
library(deSolve)

ui <- fluidPage(
  titlePanel("Representative Agent Model"),
  
  h4("Optimal Control Problem with Capital and Consumption Dynamics"),
  
  withMathJax(
    p("Maximize: \\( \\int_0^{\\infty} e^{-(\\rho - n) t} u(c(t)) dt \\)"),
    p("Subject to: \\( k'(t) = k(t)^{\\alpha} - c(t) - (\\delta + n) k(t), \\quad k(0) = k_0 \\)"),
    p("With \\( u(c) = \\frac{c^{1 - \\theta}}{1 - \\theta} \\)"),
    p("Capital Dynamics Equation: \\( k'(t) = k(t)^{\\alpha} - c(t) - (\\delta + n) k(t) \\)"),
    p("Consumption Dynamics Equation: \\( c'(t) = \\frac{c(t)}{-\\theta} (\\rho + \\delta - \\alpha k(t)^{\\alpha - 1}) \\)")
  ),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("rho", "Discount Rate (ρ)", min = 0, max = 1, value = 0.04, step = 0.01),
      sliderInput("n", "Population Growth Rate (n)", min = 0, max = 1, value = 0.01, step = 0.01),
      sliderInput("alpha", "Capital Elasticity (α)", min = 0, max = 1, value = 0.33, step = 0.01),
      sliderInput("delta", "Depreciation (δ)", min = 0, max = 1, value = 0.05, step = 0.01),
      sliderInput("theta", "Risk Aversion (θ)", min = 0.01, max = 5, value = 0.5, step = 0.01),
      sliderInput("ymax", "Maximum Consumption in Plot", min = 0, max = 10, value = 2, step = 0.1),
      sliderInput("kmax", "Maximum Capital in Plot", min = 10, max = 100, value = 70, step = 1),
      textInput("k0", "Initial Capital Value (k(0))", value = "70"),
      textInput("c0", "Initial Consumption Value (c(0))", value = "1"),
      checkboxInput("showVectors", "Show Vector Fields", value = FALSE),
      checkboxInput("showTrajectory", "Show Trajectory", value = FALSE)  # Default set to FALSE
    ),
    
    mainPanel(
      plotOutput("phasePlot"),
      h4("Steady State:"),
      verbatimTextOutput("steadyState"),
      h4("Stability Classification"),
      verbatimTextOutput("stability")
    )
  )
)

server <- function(input, output) {
  
  output$phasePlot <- renderPlot({
    rho <- input$rho
    n <- input$n
    alpha <- input$alpha
    delta <- input$delta
    theta <- input$theta
    ymax <- input$ymax
    kmax <- input$kmax  # Now kmax is controllable via the slider
    
    k0 <- as.numeric(input$k0)
    c0 <- as.numeric(input$c0)
    
    k_star <- ((delta + rho) / alpha)^(1 / (alpha - 1))
    c_star <- k_star^alpha - (delta + n) * k_star
    
    k_vals <- seq(0.1, kmax, length.out = 500)
    c_kdot0 <- k_vals^alpha - (delta + n) * k_vals
    df_kdot0 <- data.frame(k = c(0, k_vals), c = c(0, c_kdot0), curve = "k' = 0")
    df_cdot0 <- data.frame(k = rep(k_star, 500), c = seq(0, ymax, length.out = 500), curve = "c' = 0")
    df_total <- rbind(df_kdot0, df_cdot0)
    
    p <- ggplot(df_total, aes(x = k, y = c)) +
      geom_line(aes(color = curve), size = 1.5) +
      geom_point(aes(x = k_star, y = c_star), color = "black", size = 4) +
      labs(
        title = "Phase Diagram",
        x = "Capital per capita (k)",
        y = "Consumption per capita (c)",
        color = "Equilibrium Curves"
      ) +
      scale_color_manual(values = c("blue", "red"), labels = c("c' = 0", "k' = 0")) +
      scale_y_continuous(limits = c(0, ymax)) +
      theme_minimal(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    
    # Add trajectory if needed
    if (input$showTrajectory) {
      # Define the system of differential equations
      system <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          k_dot <- k^alpha - c - (delta + n) * k
          c_dot <- (c / -theta) * (rho + delta - alpha * k^(alpha - 1))
          list(c(k_dot, c_dot))
        })
      }
      
      # Initial conditions
      state <- c(k = k0, c = c0)
      parameters <- c(rho = rho, n = n, alpha = alpha, delta = delta, theta = theta)
      
      # Time points for the solution
      times <- seq(0, 100, by = 0.1)
      
      # Solve the system
      out <- ode(y = state, times = times, func = system, parms = parameters)
      trajectory_data <- as.data.frame(out)
      
      p <- p + 
        geom_path(data = trajectory_data, aes(x = k, y = c), color = "gray", linetype = "dashed") + 
        geom_point(aes(x = k0, y = c0), color = "gray", size = 4)
    }
    
    # Add vector fields if needed
    if (input$showVectors) {
      k_grid <- seq(0.1, kmax, length.out = 20)
      c_grid <- seq(0.1, ymax, length.out = 20)
      grid <- expand.grid(k = k_grid, c = c_grid)
      
      grid$k_dot <- with(grid, k^alpha - c - (delta + n) * k)
      grid$c_dot <- with(grid, c / -theta * (rho + delta - alpha * k^(alpha - 1)))
      
      mag <- sqrt(grid$k_dot^2 + grid$c_dot^2)
      grid$k_dot <- grid$k_dot / mag * 0.2
      grid$c_dot <- grid$c_dot / mag * 0.2
      
      p <- p + geom_segment(data = grid,
                            aes(x = k, y = c, xend = k + k_dot, yend = c + c_dot),
                            arrow = arrow(length = unit(0.15, "cm")),
                            color = "blue", alpha = 0.6)
    }
    
    print(p)
  })
  
  output$steadyState <- renderPrint({
    rho <- input$rho
    n <- input$n
    alpha <- input$alpha
    delta <- input$delta
    
    k_star <- ((delta + rho) / alpha)^(1 / (alpha - 1))
    c_star <- k_star^alpha - (delta + n) * k_star
    
    cat(sprintf("k* = %.4f\nc* = %.4f", k_star, c_star))
  })
  
  output$stability <- renderPrint({
    rho <- input$rho
    n <- input$n
    alpha <- input$alpha
    delta <- input$delta
    theta <- input$theta
    
    k_star <- ((delta + rho) / alpha)^(1 / (alpha - 1))
    c_star <- k_star^alpha - (delta + n) * k_star
    
    # Partial derivatives
    a11 <- (rho + delta - alpha * k_star^(alpha - 1)) / (-theta)
    a12 <- c_star * alpha * (alpha - 1) * k_star^(alpha - 2) / theta
    a21 <- -1
    a22 <- alpha * k_star^(alpha - 1) - (delta + n)
    
    # Jacobian matrix with new labels
    J <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
    colnames(J) <- c("∂c'", "∂k'")  # Column labels
    rownames(J) <- c("∂c", "∂k")  # Row labels
    
    det_J <- det(J)
    trace_J <- sum(diag(J))
    
    print(J)
    cat(sprintf("\nDeterminant: %.4f", det_J))
    cat(sprintf("\nTrace: %.4f\n", trace_J))
    
    if (det_J < 0) {
      cat("Classification: Saddle Point")
    } else if (det_J > 0 && trace_J < 0) {
      cat("Classification: Stable")
    } else if (det_J > 0 && trace_J > 0) {
      cat("Classification: Unstable")
    } else {
      cat("Classification: Indeterminate")
    }
  })
}

shinyApp(ui, server)
