maxLikOverlap.PH<- 
  function (ellipse1, ellipse2, siber.object, p.interval = 0.95, 
          n = 100, do.plot = FALSE) 
{
  tmp <- ellipse1
  c.1 <- 'AllTaxa'
  e.1 <- tmp
  coords.1 <- addEllipse(siber.object$ML.mu[[c.1]][, , e.1], 
                         siber.object$ML.cov[[c.1]][, , e.1], m = siber.object$sample.sizes[c.1, 
                                                                                            e.1], small.sample = TRUE, n = n, p.interval = p.interval, 
                         ci.mean = FALSE, do.plot = FALSE)
  area.1 <- siberConvexhull(coords.1[, 1], coords.1[, 2])$TA
  tmp <- ellipse2
  c.2 <- 'AllTaxa'
  e.2 <- tmp
  coords.2 <- addEllipse(siber.object$ML.mu[[c.2]][, , e.2], 
                         siber.object$ML.cov[[c.2]][, , e.2], m = siber.object$sample.sizes[c.2, 
                                                                                            e.2], small.sample = TRUE, n = n, p.interval = p.interval, 
                         ci.mean = FALSE, do.plot = FALSE)
  area.2 <- siberConvexhull(coords.2[, 1], coords.2[, 2])$TA
  overlap <- abs(spatstat.utils::overlap.xypolygon(list(x = coords.1[, 
                                                                     1], y = coords.1[, 2]), list(x = coords.2[, 1], y = coords.2[, 
                                                                                                                                  2])))
  out <- c(area.1 = area.1, area.2 = area.2, overlap = overlap)
  return(out)
  }
