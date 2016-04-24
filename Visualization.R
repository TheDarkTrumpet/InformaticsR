library(ggplot2)

# mtcars = test data set.
qplot(mtcars$mpg)

qplot(cyl, mpg, data=mtcars)
qplot(mpg, qsec, data=mtcars)

qplot(mpg, data=mtcars, geom="density")     # Line Graph
qplot(wt, mpg, data=mtcars, geom="smooth")  # Smooth Line (Ribbon)

qplot(wt, mpg, data=mtcars, geom=c("smooth", "point"))  # Smooth Line (Ribbon)
qplot(mpg, qsec, data=mtcars, geom=c("smooth", "point"))  # Smooth Line (Ribbon)

qplot(wt, mpg, data=mtcars, geom=c("boxplot", "jitter"), fill=gear)

# Make a copy of mtcars
mtCars2 <- mtcars
mtCars2$cyl <- factor(mtcars$cyl, levels=c(4,6,8), labels=c("4cyl", "6cyl", "8cyl"))

mtcars$cyl
mtCars2$cyl  # Replaced 4 = 4cyl, 6 = 6cyl, 8 = 8cyl -- did this above

# Various of the same general plot
qplot(cyl, mpg, data=mtCars2, geom="boxplot")
qplot(cyl, mpg, data=mtCars2, geom=c("boxplot", "jitter"), fill=gear)
qplot(cyl, mpg, data=mtCars2, geom=c("boxplot", "jitter"), fill=gear)
qplot(cyl, mpg, data=mtCars2, geom=c("boxplot", "jitter"), fill=cyl)
