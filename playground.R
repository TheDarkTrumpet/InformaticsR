# Some of this information is taken from RSeminarScript, and some from Visualization
# it is more playing around with data

df1 = data.frame(id=1:6,
  gender=factor(rep(c("M","F"),3)),
  treatment1=factor(rep(LETTERS[1:3],2)),
  treatment2=factor(rep(LETTERS[1:2],each=3)),
  response=rnorm(6))

# Selecting information from the dataframe (All do the same thing)
df1[df1$treatment1 == 'C',]
df1[c(3,6),]
df1[with(df1, which(treatment1=="C")),]