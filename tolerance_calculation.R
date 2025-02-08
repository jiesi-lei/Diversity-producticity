# calculate number of tolerant species for each combination
comm = read.csv("../raw_data/occurrence.csv",row.names = 1)
toler =  read.csv("data/rrn_tolerance/rrn.csv",row.names = 1)

salt.tol <- starve.tol <- c()
for (i in 1:81){
salt.tol = c(salt.tol, sum(comm[i,]*toler$Salinity.tolerant))
starve.tol = c(starve.tol, sum(comm[i,]*toler$Starvation.tolerant))
}

copy_number = data.frame(rownames(comm),salt.tol,starve.tol)
write.csv(copy_number, "data/rrn_tolerance/tolerance.csv")
