plotNowcast <- function(data, nowcast, title = "Nowcast") {
  #
  # Descrition
  # Plot nowcast
  #
  # Arguments
  # data     Dataframe, output from dataSetup
  # nowcast  List, output from nowcast
  # title    Title of the plot
  #
  # Value
  # ggplot object
  
  # Prepare data for nowcast plot
  tmp1 <- data %>%
    mutate(Reported = Reported %>%
        recode_factor("Future" = "Not yet reported") %>% 
        factor(
          levels = c("Reported", "Not yet reported", "Nowcast"),
          labels = c("Reported", "Not yet reported", "Nowcast"))) %>% 
    group_by(Date, Reported) %>%
    summarize(Cases = sum(Cases))
  tmp2 <- nowcast.list$nowcast %>% 
    filter(
      Date %in% ((
        data %>%
          filter(Reported == "Not yet reported"))$Date %>%
          unique %>%
          sort))
  
  # Make the plot
  plot <- ggplot(
    data = tmp1,
    mapping = aes(x = Date, y = Cases, fill = Reported)) +
    geom_col(
      width = 1,
      position = position_stack(reverse = TRUE)) +
    geom_step(
      data = tmp2,
      mapping = aes(x = Date + 0.5, y = med),
      colour = oran,
      size = 0.5,
      direction = "vh",
      inherit.aes = FALSE) +
    geom_crossbar(
      data = tmp2,
      mapping = aes(x = Date, y = med, ymin = lwr, ymax = upr),
      fill = adjustcolor(oran, alpha = 0.3), colour = NA, width = 1,
      inherit.aes = FALSE) +
    scale_fill_manual(
      values = c(blue, grey, oran),
      name = "",
      drop = FALSE) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b %d",
      expand = c(0, 0)) +
    scale_y_continuous(
      limits = c(0, with(data, max(tapply(Cases, Date, sum)) + 1)),
      expand = c(0, 0),
      oob = squish) +
    labs(
      x = "Time of symptoms onset",
      y = "Number of symptomatic cases",
      title = paste(title, "-", format(with(subset(data, Reported == "Reported"), max(Date)), format = "%b %d, %Y"))) +
    theme_bw() +
    theme(
      plot.margin = margin(t = 0.1, r = 0.5, b = 0.1, l = 0.1, unit = "cm"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top")
  
  # Return output
  return(plot)
  
}
