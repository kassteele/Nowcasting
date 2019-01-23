plotDelayDist <- function(data, nowcast, title = "Reporting delay distribution") {
  #
  # Description
  # Plot reporting delay distribution
  #
  # Arguments
  # data     Dataframe, output from dataSetup
  # nowcast  List, output from nowcast
  # title    Title of the plot
  #
  # Value
  # ggplot object
  
  # Prepare data for delay distribution plot
  tmp <- right_join(
    x = nowcast$f.delay,
    y = data %>% select(Date, Delay, Reported))

  # Prepare contour labels data for delay distribution plot
  contourlabels.data <- lapply(
    X = with(nowcast$f.delay, contourLines(
      x = unique(Date),
      y = unique(Delay),
      z = matrix(f.delay, length(unique(Date)), length(unique(Delay))),
      levels = seq(from = 0, to = 0.20, by = 0.01))),
    FUN = function(l) {
      with(l, {
        i <- which.max(x)
        return(data.frame(level = level, x = as.Date(x[i], origin = "1970-01-01"), y = y[i]))
      })
    }) %>% 
    do.call(what = "rbind")
  
  # Delay distribution plot
  plot <- ggplot(
    data = tmp,
    mapping = aes(x = Date, y = Delay, fill = rescale(f.delay) + 100*(as.numeric(Reported) - 1), z = f.delay)) +
    geom_raster() +
    geom_contour(
      binwidth = 0.01,
      size = 0.2,
      colour = "black") +
    geom_text(
      data = contourlabels.data,
      mapping = aes(x = x, y = y, label = level),
      size = 3,
      nudge_x = 5,
      check_overlap = TRUE,
      inherit.aes = FALSE) +
    scale_fill_gradientn(
      limits = c(0, 101),
      colours = c(blue.pal(n = 5), oran.pal(n = 5)),
      values = rescale(c(outer(seq(from = 0, to = 1, length = 5), c(0, 100), "+"))),
      na.value = NA) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b %d",
      expand = c(0, 0)) +
    scale_y_continuous(
      expand = c(0, 0)) +
    labs(
      x = "Time of symptoms onset",
      y = "Reporting delay",
      title = paste(title, "-", format(max(nowcast$f.delay$Date), "%b %d, %Y"))) +
    coord_fixed() +
    theme_bw() +
    theme(
      plot.margin = margin(t = 0.1, r = 0.5, b = 0.1, l = 0.1, unit = "cm"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none")
  
  # Return outpout
  return(plot)
  
}
