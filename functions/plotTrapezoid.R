plotTrapezoid <- function(data, title = "Reporting trapezoid") {
  #
  # Descrition
  # Plot reporting trapezoid
  #
  # Arguments
  # data   Dataframe, output from dataSetupt
  # title  Title of the plot
  #
  # Value
  # ggplot object

  # Make the plot
  plot <- ggplot(
    data = data,
    mapping = aes(x = Date, y = Delay, fill = rescale(log1p(Cases)) + 100*(as.numeric(Reported) - 1))) +
    geom_raster() +
    scale_fill_gradientn(
      limits = c(0, 201),
      colours = c(blue.pal(n = 5), oran.pal(n = 5), grey.pal(n = 5)),
      values = rescale(c(outer(seq(from = 0, to = 1, length = 5), c(0, 100, 200), "+")))) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b %d",
      expand = c(0, 0)) +
    scale_y_continuous(
      expand = c(0, 0)) +
    coord_fixed(
      ratio = 1) +
    labs(
      x = "Time of symptoms onset",
      y = "Reporting delay",
      title = paste(title, "-", format(with(subset(data, Reported == "Reported"), max(Date)), format = "%b %d, %Y"))) +
    theme_bw() +
    theme(
      plot.margin = margin(t = 0.1, r = 0.5, b = 0.1, l = 0.1, unit = "cm"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none")
  
  # Return output
  return(plot)
  
}
