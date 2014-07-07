tf <- data.frame(temperature = seq(15, 31, 0.01), prob = dnorm(seq(15, 31, 0.01), 23, 2))
p <- ggplot(tf, aes(x = temperature, y = prob))
p <- p + geom_line(lwd = 2)
p <- p + theme_bw(20)
p <- p + scale_x_continuous("Temperature tomorrow")
p <- p + scale_y_continuous("Probability density")
p
ggsave("weather.pdf", p)
