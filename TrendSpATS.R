

Trend.SpATS <-
  function(x, main = NULL, depict.missing = FALSE, graphics=F  ,...) {
    
    all.in.one = TRUE
    annotated = FALSE
    
    xlab <- x$terms$spatial$terms.formula$x.coord
    ylab <- x$terms$spatial$terms.formula$y.coord
    x.coord <- x$data[,xlab]
    y.coord <- x$data[,ylab]
    
    if(is.null(main)) main = paste("Trait: ", x$model$response, sep = "")
    
    columns <- seq(min(x.coord), max(x.coord), by = min(diff(sort(unique(x.coord)))))
    rows <- seq(min(y.coord), max(y.coord), by = min(diff(sort(unique(y.coord)))))
    
    setNumericRounding(2)
    
    xy.coord <- data.table(expand.grid(columns = columns, rows = rows))
    setkeyv(xy.coord, c("rows", "columns"))
    ONE <- rep(1, length(x.coord))    
    df <- data.table(columns = x.coord, rows = y.coord, weights = x$data$weights, ONE = ONE)
    setkeyv(df, c("rows", "columns"))
    df <- df[xy.coord]
    df <- df[order(df$columns, df$rows),]
    if(depict.missing) {
      df$ONE[is.na(df$ONE)] <- 1
      }  else {
      df$ONE[df$weights == 0 | is.na(df$weights)] <- NA 
      }
    
    p1 <- if(length(columns) > 100) 1 else 100%/%length(columns) + 1
    p2 <- if(length(rows) > 100) 1 else 100%/%length(rows) + 1
    
    fit.spatial.trend <- obtain.spatialtrend(x, grid = c(length(columns)*p1, length(rows)*p2))
    Mf = kronecker(matrix(df$ONE, ncol = length(columns), nrow = length(rows)), matrix(1, p2, p1))
    
    colors = topo.colors(100)
    
    
    if (graphics) {
      main.legends <- c('Raw data', 'Fitted data', 'Residuals', 'Fitted Spatial Trend', ifelse(x$model$geno$as.random, "Genotypic BLUPs", "Genotypic BLUEs"), 'Histogram')
      if(all.in.one) {
        op <- par( oma = c(ifelse(annotated, 12, 2), 1, 3, 2), mar = c(2.7, 4, 2.5, 2.5), mgp = c(1.7, 0.5, 0))                
      } else {
        if(!is.null(main))
          main.legends <- rep(main, length(main.legends))
      }
      
      fields::image.plot(fit.spatial.trend$col.p, fit.spatial.trend$row.p, t(fit.spatial.trend$fit*Mf), main = main.legends[4], col = colors, xlab = xlab, ylab = ylab, graphics.reset = F, ...)
      
      if(all.in.one) {
        title("")
        mtext(main, cex = 1.5, outer = TRUE, side = 3)
        par(op)
      }
    }
    
    Spatial <- list(Col.p=fit.spatial.trend$col.p, Row.p=fit.spatial.trend$row.p, GL=t(fit.spatial.trend$fit*Mf) )
    invisible(Spatial)
    
  }

# Example
# Trend.SpATS(Model,depict.missing = F,graphics = F)GL


